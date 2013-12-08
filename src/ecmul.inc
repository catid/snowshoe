// Elliptic curve point multiplication

#include "ecpt.inc"
#include "misc.inc"
#include "recode.inc"

/*
 * Multiplication by generator point using LSB-set comb method [1] with w=6,v=7
 *
 * Preconditions:
 * 	0 < k < q
 *
 * Multiplies the point by k and stores the result in R
 */

// R = k * G
static void ec_mul_gen(const u64 k[4], ecpt_affine &R) {
	// Recode scalar
	u64 rk[4];
	ec_recode_comb(k, rk);

	// Unroll first evaluation loop
	ecpt T[MG_n];
	ec_table_select_comb_gen(rk, MG_s - 1, T);
	fe_set_smallk(1, T[0].z);

	// X = T[0] + T[1] + ...
	ufe t2b;
	ecpt X;
	ec_add(T[0], T[1], X, true, true, false, t2b);
	for (int jj = 2; jj < MG_n; ++jj) {
		ec_add(X, T[jj], X, true, false, false, t2b);
	}

	// Evaluate
	for (int ii = MG_s - 2; ii >= 0; --ii) {
		ec_table_select_comb_gen(rk, ii, T);

		ec_dbl(X, X, false, t2b);
		for (int jj = 0; jj < MG_n; ++jj) {
			ec_add(X, T[jj], X, true, false, false, t2b);
		}
	}

	// Compute affine coordinates in R
	ec_affine(X, R);
}

/*
 * Multiplication by variable base point using GLV-SAC method [1] with m=2
 *
 * Preconditions:
 * 	0 < k < q
 *
 * Multiplies the point by k * 4 and stores the result in R
 */

// R = k*4*P
static void ec_mul(const u64 k[4], const ecpt_affine &P0, ecpt_affine &R) {
	// Decompose scalar into subscalars
	ufp a, b;
	s32 asign, bsign;
	gls_decompose(k, asign, a, bsign, b);

	// Q0 = endomorphism of P0
	ecpt_affine Q0;
	gls_morph(P0.x, P0.y, Q0.x, Q0.y);

	// Set base point sign
	ec_cond_neg_affine(bsign, Q0);

	// Expand P, Q to extended coordinates
	ecpt P, Q;
	ec_expand(P0, P);
	ec_expand(Q0, Q);

	// Set base point sign
	ec_cond_neg(asign, P);

	// Precompute multiplication table
	ecpt table[8];
	ec_gen_table_2(P, Q, table);

	// Recode subscalars
	u32 recode_bit = ec_recode_scalars_2(a, b, 128);

	// Initialize working point
	ecpt X;
	ec_table_select_2(table, a, b, 126, true, X);

	// Evaluate
	ufe t2b;
	for (int ii = 124; ii >= 0; ii -= 2) {
		ecpt T;
		ec_table_select_2(table, a, b, ii, true, T);

		ec_dbl(X, X, false, t2b);
		ec_dbl(X, X, false, t2b);
		ec_add(X, T, X, false, false, false, t2b);
	}

	// If bit == 1, X <- X + P (inverted logic from [1])
	ec_cond_add(recode_bit, X, P, X, true, false, t2b);

	// Multiply by 4 to avoid small subgroup attack
	ec_dbl(X, X, false, t2b);
	ec_dbl(X, X, false, t2b);

	// Compute affine coordinates in R
	ec_affine(X, R);
}

/*
 * Simultaneous multiplication by two base points,
 * where one is variable and the other is the generator point,
 * using a combination of GLV-SAC with m=2 and LSB-set comb with w=8,v=1 [1].
 *
 * Preconditions:
 * 	0 < a,b < q
 *
 * Multiplies the result of aG + bP by 4 and stores it in R
 */

// R = a*4*G + b*4*P
static void ec_simul_gen(const u64 a[4], const u64 b[4], const ecpt_affine &P0, ecpt_affine &R) {
	// Decompose scalar into subscalars
	ufp b1, b2;
	s32 b1sign, b2sign;
	gls_decompose(b, b1sign, b1, b2sign, b2);

	// Q0 = endomorphism of P0
	ecpt_affine Q0;
	gls_morph(P0.x, P0.y, Q0.x, Q0.y);

	// Set base point sign
	ec_cond_neg_affine(b2sign, Q0);

	// Expand base points to extended coordinates
	ecpt P, Q;
	ec_expand(P0, P);
	ec_expand(Q0, Q);

	// Set base point sign
	ec_cond_neg(b1sign, P);

	// Precompute multiplication table
	ecpt qtable[8];
	ec_gen_table_2(P, Q, qtable);

	// Recode subscalars
	u64 a1[4];
	ec_recode_comb(a, a1);
	u32 recode_bit = ec_recode_scalars_2(b1, b2, 128);

	// Initialize working point
	ecpt X;
	ec_table_select_2(qtable, b1, b2, 126, false, X);

	// Evaluate
	ufe t2b;
	ecpt T;
	for (int ii = 124; ii >= 32; ii -= 2) {
		ec_table_select_2(qtable, b1, b2, ii, false, T);

		ec_dbl(X, X, false, t2b);
		ec_dbl(X, X, false, t2b);
		ec_add(X, T, X, false, false, false, t2b);
	}

	// For the last 32 doubles, interleave ec_mul_gen adds
	for (int ii = 30; ii >= 0; ii -= 2) {
		ec_dbl(X, X, false, t2b);

		ec_table_select_comb_81(a1, ii+1, T);
		ec_add(X, T, X, true, false, false, t2b);

		ec_dbl(X, X, false, t2b);

		ec_table_select_comb_81(a1, ii, T);
		ec_add(X, T, X, true, false, false, t2b);

		ec_table_select_2(qtable, b1, b2, ii, false, T);
		ec_add(X, T, X, false, false, false, t2b);
	}

	// If bit == 1, X <- X + P1 (inverted logic from [1])
	if (recode_bit != 0) {
		ec_add(X, P, X, true, false, false, t2b);
	}

	// Multiply by 4 to avoid small subgroup attack
	ec_dbl(X, X, false, t2b);
	ec_dbl(X, X, false, t2b);

	// Compute affine coordinates in R
	ec_affine(X, R);
}

/*
 * Simultaneous multiplication by two variable base points
 * using GLV-SAC with m=4 [1].
 *
 * Preconditions:
 * 	0 < a,b < q
 *
 * Multiplies the result of aP + bQ by 4 and stores it in R
 */

// R = a*4*P + b*4*Q
static void ec_simul(const u64 a[4], const ecpt_affine &P, const u64 b[4], const ecpt_affine &Q, ecpt_affine &R) {
	// Decompose scalar into subscalars
	ufp a0, a1, b0, b1;
	s32 a0sign, a1sign, b0sign, b1sign;
	gls_decompose(a, a0sign, a0, a1sign, a1);
	gls_decompose(b, b0sign, b0, b1sign, b1);

	// Compute endomorphism of base points
	ecpt_affine Pe, Qe;
	gls_morph(P.x, P.y, Pe.x, Pe.y);
	gls_morph(Q.x, Q.y, Qe.x, Qe.y);

	// Set base point signs
	ec_cond_neg_affine(a1sign, Pe);
	ec_cond_neg_affine(b1sign, Qe);

	// Expand base points
	ecpt P0, P1, Q0, Q1;
	ec_expand(P, P0);
	ec_expand(Pe, P1);
	ec_expand(Q, Q0);
	ec_expand(Qe, Q1);

	// Set base point signs
	ec_cond_neg(a0sign, P0);
	ec_cond_neg(b0sign, Q0);

	// Precompute multiplication table
	ecpt table[8];
	ec_gen_table_4(P0, P1, Q0, Q1, table);

	// Recode scalar
	u32 recode_bit = ec_recode_scalars_4(a0, a1, b0, b1, 127);

	// Initialize working point
	ecpt X;
	ec_table_select_4(table, a0, a1, b0, b1, 126, X);

	// Evaluate
	ufe t2b;
	for (int ii = 125; ii >= 0; --ii) {
		ecpt T;
		ec_table_select_4(table, a0, a1, b0, b1, ii, T);

		ec_dbl(X, X, false, t2b);
		ec_add(X, T, X, false, false, false, t2b);
	}

	// If bit == 1, X <- X + P (inverted logic from [1])
	ec_cond_add(recode_bit, X, P0, X, true, false, t2b);

	// Multiply by 4 to avoid small subgroup attack
	ec_dbl(X, X, false, t2b);
	ec_dbl(X, X, false, t2b);

	// Compute affine coordinates in R
	ec_affine(X, R);
}

