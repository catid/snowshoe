// Elliptic curve point multiplication

#include "ecpt.inc"
#include "misc.inc"
#include "recode.inc"

/*
 * Multiplication by generator point using LSB-set comb method [1] with w=6,v=7
 *
 * Note that this function actually does support k=0, and it will return
 * the same as k=q in this case; the result is the identity element (0, 1).
 *
 * Optionally adds point 4E.  Pass null to skip that step.
 *
 * Preconditions:
 * 	0 < k < q
 *
 * Multiplies the point by k and stores the result in R
 */

// R = kG
static void ec_mul_gen(const u64 k[4], ecpt &R, ufe &r2b) {
	// Recode scalar
	u64 kp[4];
	u32 recode_lsb = ec_recode_scalar_comb_gen(k, kp);

	// Unroll first evaluation loop
	ecpt T[MG_v];
	ec_table_select_comb_gen(kp, MG_e - 1, T);
	fe_set_smallk(1, T[0].z);

	// X = T[0] + T[1] + T[2]
	ufe t2b;
	ecpt X;
	ec_add(T[0], T[1], X, true, true, false, t2b);
	for (int jj = 2; jj < MG_v; ++jj) {
		ec_add(X, T[jj], X, true, false, false, t2b);
	}

	// Evaluate
	for (int ii = MG_e - 2; ii >= 0; --ii) {
		ec_table_select_comb_gen(kp, ii, T);

		ec_dbl(X, X, false, t2b);
		for (int jj = 0; jj < MG_v; ++jj) {
			ec_add(X, T[jj], X, true, false, false, t2b);
		}
	}

	// NOTE: Do conditional addition here rather than after the ec_cond_neg
	// (this is an error in the paper)
	// If carry bit is set, add 2^(w*d)
	ec_cond_add((kp[3] >> 60) & 1, X, *GEN_FIX, X, true, false, t2b);

	// If recode_lsb == 1, R = -R
	ec_cond_neg(recode_lsb, X, R);

	// Copy t2b out
	fe_set(t2b, r2b);
}

/*
 * Multiplication by variable base point using GLV-SAC method [1] with m=2
 *
 * Note that this function will fail if k=0.
 *
 * Preconditions:
 * 	0 < k < q
 *
 * Multiplies the point by k and stores the result in R, r2b
 */

// R = kP
static void ec_mul(const u64 k[4], const ecpt &P0, bool z1, ecpt &R, ufe &r2b) {
	// Decompose scalar into subscalars
	ufp a, b;
	s32 asign, bsign;
	gls_decompose(k, asign, a, bsign, b);

	// Q = endomorphism(P)
	ecpt P, Q;
	if (z1) {
		// Simplify the endomorphism calculation
		gls_morph(P0.x, P0.y, Q.x, Q.y);
		fe_mul(Q.x, Q.y, Q.t);
		fe_set_smallk(1, Q.z);
	} else {
		gls_morph_ext(P0, Q);
	}

	// Set base point signs
	ec_cond_neg(asign, P0, P);
	ec_cond_neg(bsign, Q, Q);

	// Precompute multiplication table
	ecpt table[8];
	ec_gen_table_2(P, Q, z1, table);

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

	// If bit == 1, R <- X + P (inverted logic from [1])
	ec_cond_add(recode_bit, X, P, R, z1, false, t2b);

	// Copy t2b out
	fe_set(t2b, r2b);
}

/*
 * Simultaneous multiplication by two base points,
 * where one is variable and the other is the generator point,
 * using a combination of GLV-SAC with m=2 and LSB-set comb with w=8,v=1 [1].
 *
 * Note that this function will fail if b=0.
 * Note that this function actually does support a=0, and it will return
 * the same as a=q in this case; the result = bP.
 *
 * Preconditions:
 * 	0 < a,b < q
 *
 * Performs aG + bP and stores it in R, r2b
 */

// R = aG + bP
static void ec_simul_gen(const u64 a[4], const u64 b[4], const ecpt &P0, bool z1, ecpt &R, ufe &r2b) {
	// Decompose scalar into subscalars
	ufp b1, b2;
	s32 b1sign, b2sign;
	gls_decompose(b, b1sign, b1, b2sign, b2);

	// Q = endomorphism(P)
	ecpt P, Q;
	if (z1) {
		// Simplify the endomorphism calculation
		gls_morph(P0.x, P0.y, Q.x, Q.y);
		fe_mul(Q.x, Q.y, Q.t);
		fe_set_smallk(1, Q.z);
	} else {
		gls_morph_ext(P0, Q);
	}

	// Set base point signs
	ec_cond_neg(b1sign, P0, P);
	ec_cond_neg(b2sign, Q, Q);

	// Precompute multiplication table
	ecpt qtable[8];
	ec_gen_table_2(P, Q, z1, qtable);

	// Recode subscalars
	u64 a1[4];
	const u32 comb_lsb = ec_recode_scalar_comb_81(a, a1);
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

		ec_table_select_comb_81(comb_lsb, a1, ii+1, T);
		ec_add(X, T, X, true, false, false, t2b);

		ec_dbl(X, X, false, t2b);

		ec_table_select_comb_81(comb_lsb, a1, ii, T);
		ec_add(X, T, X, true, false, false, t2b);

		ec_table_select_2(qtable, b1, b2, ii, false, T);
		ec_add(X, T, X, false, false, false, t2b);
	}

	// If bit == 1, R <- X + P (inverted logic from [1])
	ec_cond_add(recode_bit, X, P, R, z1, false, t2b);

	// Copy t2b out
	fe_set(t2b, r2b);
}

/*
 * Simultaneous multiplication by two variable base points
 * using GLV-SAC with m=4 [1].
 *
 * Note that this function will fail if a,b=0.
 *
 * Preconditions:
 * 	0 < a,b < q
 *
 * Performs aP + bQ and stores the result in R
 */

// R = aP + bQ
static void ec_simul(const u64 a[4], const ecpt &P0, bool pz1, const u64 b[4], const ecpt &Q0, bool qz1, ecpt &R, ufe &r2b) {
	// Decompose scalar into subscalars
	ufp a0, a1, b0, b1;
	s32 a0sign, a1sign, b0sign, b1sign;
	gls_decompose(a, a0sign, a0, a1sign, a1);
	gls_decompose(b, b0sign, b0, b1sign, b1);

	// Q = endomorphism(P)
	ecpt P, Pe, Q, Qe;
	if (pz1) {
		// Simplify the endomorphism calculation
		gls_morph(P0.x, P0.y, Pe.x, Pe.y);
		fe_mul(Pe.x, Pe.y, Pe.t);
		fe_set_smallk(1, Pe.z);
	} else {
		gls_morph_ext(P0, Pe);
	}
	if (qz1) {
		// Simplify the endomorphism calculation
		gls_morph(Q0.x, Q0.y, Qe.x, Qe.y);
		fe_mul(Qe.x, Qe.y, Qe.t);
		fe_set_smallk(1, Qe.z);
	} else {
		gls_morph_ext(Q0, Qe);
	}

	// Set base point signs
	ec_cond_neg(a0sign, P0, P);
	ec_cond_neg(a1sign, Pe, Pe);
	ec_cond_neg(b0sign, Q0, Q);
	ec_cond_neg(b1sign, Qe, Qe);

	// Precompute multiplication table
	ecpt table[8];
	ec_gen_table_4(P, Pe, pz1, Q, Qe, qz1, table);

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

	// If bit == 1, R <- X + P (inverted logic from [1])
	ec_cond_add(recode_bit, X, P, R, pz1, false, t2b);

	// Copy t2b out
	fe_set(t2b, r2b);
}

