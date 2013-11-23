// Elliptic curve point operations
#include "ecpt.cpp"

#include <iostream>
#include <iomanip>
using namespace std;

#include "Clock.hpp"

/*
 * Mask a random number to produce a compatible scalar for multiplication
 */

void ec_mask_scalar(u64 k[4]) {
	// Prime order of the curve = q, word-mapped:
	// 0x0FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFA6261414C0DC87D3CE9B68E3B09E01A5
	//   (      3       )(       2      )(       1      )(       0      )

	// Clear high 5 bits
	// Clears one extra bit to simplify key generation
	k[3] &= 0x07FFFFFFFFFFFFFFULL;

	// Largest value after filtering:
	// 0x07FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
	//   (      3       )(       2      )(       1      )(       0      )
}

/*
 * GLV-SAC Scalar Recoding Algorithm for m=2 [1]
 *
 * Returns low bit of 'a'.
 */

static u32 ec_recode_scalars_2(ufp &a, ufp &b, const int len) {
	u32 lsb = ((u32)a.i[0] & 1) ^ 1;
	a.w -= lsb;
	a.w >>= 1;
	a.w |= (u128)1 << (len - 1);

	const u128 an = ~a.w;
	u128 mask = 1;
	for (int ii = 1; ii < len; ++ii) {
		const u128 anmask = an & mask;
		b.w += (b.w & anmask) << 1;
		mask <<= 1;
	}

	return lsb;
}

/*
 * GLV-SAC Scalar Recoding Algorithm for m=4 [1]
 *
 * Returns low bit of 'a'.
 */

static u32 ec_recode_scalars_4(ufp &a, ufp &b, ufp &c, ufp &d, const int len) {
	u32 lsb = ((u32)a.i[0] & 1) ^ 1;
	a.w -= lsb;
	a.w >>= 1;
	a.w |= (u128)1 << (len - 1);

	const u128 an = ~a.w;
	u128 mask = 1;
	for (int ii = 1; ii < len; ++ii) {
		const u128 anmask = an & mask;
		b.w += (b.w & anmask) << 1;
		c.w += (c.w & anmask) << 1;
		d.w += (d.w & anmask) << 1;
		mask <<= 1;
	}

	return lsb;
}


//// Constant-time Generator Base Multiplication

/*
 * Precomputed table generation
 *
 * Using GLV-SAC Precomputation with m=2 [1], assuming window size of 2 bits
 *
 * Table index is simply = (a0 ^ a1) || b1 || b0
 *
 * The differences above from [1] seem to improve the efficiency of evaulation
 * and they make the code easier to analyze.
 */

void ec_gen_table_128(const ecpt &a, const ecpt &b, ecpt TABLE[128]) {
	for (u32 jj = 0; jj < 128; ++jj) {
		u32 ii = jj;

		s32 ak = 0;
		if (ii & (1 << 4)) {
			ak += 1;
		} else {
			ak -= 1;
		}
		if (ii & (1 << 5)) {
			ak += 2;
		} else {
			ak -= 2;
		}
		if (ii & (1 << 6)) {
			ak += 4;
		} else {
			ak -= 4;
		}
		if (ii & (1 << 7)) {
			ak += 8;
		} else {
			ak -= 8;
		}
		s32 bk = 0;
		if (ii & (1 << 0)) {
			if (ii & (1 << 4)) {
				bk += 1;
			} else {
				bk -= 1;
			}
		}
		if (ii & (1 << 1)) {
			if (ii & (1 << 5)) {
				bk += 2;
			} else {
				bk -= 2;
			}
		}
		if (ii & (1 << 2)) {
			if (ii & (1 << 6)) {
				bk += 4;
			} else {
				bk -= 4;
			}
		}
		if (ii & (1 << 3)) {
			if (ii & (1 << 7)) {
				bk += 8;
			} else {
				bk -= 8;
			}
		}
		cout << ii << " : [" << ak << "]a + [" << bk << "]b";

		ii |= 0x80;
		ii ^= 0x70;

		ak = 0;
		if (ii & (1 << 4)) {
			ak += 1;
		} else {
			ak -= 1;
		}
		if (ii & (1 << 5)) {
			ak += 2;
		} else {
			ak -= 2;
		}
		if (ii & (1 << 6)) {
			ak += 4;
		} else {
			ak -= 4;
		}
		if (ii & (1 << 7)) {
			ak += 8;
		} else {
			ak -= 8;
		}
		bk = 0;
		if (ii & (1 << 0)) {
			if (ii & (1 << 4)) {
				bk += 1;
			} else {
				bk -= 1;
			}
		}
		if (ii & (1 << 1)) {
			if (ii & (1 << 5)) {
				bk += 2;
			} else {
				bk -= 2;
			}
		}
		if (ii & (1 << 2)) {
			if (ii & (1 << 6)) {
				bk += 4;
			} else {
				bk -= 4;
			}
		}
		if (ii & (1 << 3)) {
			if (ii & (1 << 7)) {
				bk += 8;
			} else {
				bk -= 8;
			}
		}
		cout << " -> [" << ak << "]a + [" << bk << "]b";

		ecpt p;
		ufe t2b;

		ec_set(a, p);
		for (int kk = 1; kk < ak; ++kk) {
			ec_add(p, a, p, false, true, true, t2b);
		}

		ecpt p2;
		ec_set(b, p2);
		if (bk < 0) {
			ec_neg(p2, p2);
			bk = -bk;
			cout << "-";
		}

		for (int kk = 0; kk < bk; ++kk) {
			ec_add(p, p2, p, false, true, true, t2b);
		}

		ec_affine(p, true, p);

		ec_set(p, TABLE[jj]);

		cout << endl;
	}
}

/*
 * Table selection rule:
 *
 * k = (a2^a3) || (a1^a3) || (a0^a3) || b3 || b2 || b1 || b0
 */

static CAT_INLINE void ec_table_select_128(const ecpt *table, const ufp &a, const ufp &b, const int index, ecpt &r) {
	u32 bits = (u32)(a.w >> index) & 15;
	u32 mask = -(s32)(bits >> 3) & 7;
	u32 k = (bits & mask) << 4;
	k |= (u32)(b.w >> index) & 15;

	cout << k << endl;

	for (int ii = 0; ii < 128; ++ii) {
		// Generate a mask that is -1 if ii == index, else 0
		const u128 mask = ec_gen_mask(ii, k);

		// Add in the masked table entry
		ec_set_mask(table[ii], mask, r);
	}

	ec_cond_neg((mask & 1) ^ 1, r);
}

/*
 * Multiplication by variable base point
 *
 * Point must be in extended projective coordinates with T and Z values.
 *
 * The resulting point has undefined T and Z values, so must be expanded with
 * ec_expand() before using as input to other math functions.
 */

// R = kP
void ec_mul_gen(const u64 k[4], ecpt &R) {
	// Decompose scalar into subscalars
	ufp a, b;
	s32 asign, bsign;
	gls_decompose(k, asign, a, bsign, b);

	// Set base point signs
	ecpt P, Q;
	ec_set(EC_G, P);
	ec_set(EC_EG, Q);
	ec_cond_neg(asign, P);
	ec_cond_neg(bsign, Q);

	// Precompute multiplication table
	ecpt table[128];
	ec_gen_table_128(P, Q, table);

	// Recode subscalars
	u32 recode_bit = ec_recode_scalars_2(a, b, 128);

	// Initialize working point
	ecpt X;
	ec_table_select_128(table, a, b, 124, X);

	ufe t2b;
	for (int ii = 120; ii >= 0; ii -= 4) {
		ecpt T;
		ec_table_select_128(table, a, b, ii, T);

		ec_dbl(X, X, false, t2b);
		ec_dbl(X, X, false, t2b);
		ec_dbl(X, X, false, t2b);
		ec_dbl(X, X, false, t2b);
		ec_add(X, T, X, false, false, false, t2b);
	}

	// If bit == 1, X <- X + P (inverted logic from [1])
	ec_cond_add(recode_bit, X, P, X, false, false, t2b);

	// Compute affine coordinates in R
	ec_affine(X, false, R);
}


//// Constant-time Variable Base Multiplication

/*
 * Precomputed table generation
 *
 * Using GLV-SAC Precomputation with m=2 [1], assuming window size of 2 bits
 *
 * Window of 2 bits table selection:
 *
 * aa bb -> evaluated (unsigned table index), sign
 * 00 00    -3a + 0b (0)-
 * 00 01    -3a - 1b (1)-
 * 00 10    -3a - 2b (2)-
 * 00 11    -3a - 3b (3)-
 * 01 00    -1a + 0b (4)-
 * 01 01    -1a + 1b (5)-
 * 01 10    -1a - 2b (6)-
 * 01 11    -1a - 1b (7)-
 * 10 00    1a + 0b (4)+
 * 10 01    1a - 1b (5)+
 * 10 10    1a + 2b (6)+
 * 10 11    1a + 1b (7)+
 * 11 00    3a + 0b (0)+
 * 11 01    3a + 1b (1)+
 * 11 10    3a + 2b (2)+
 * 11 11    3a + 3b (3)+
 *
 * Table index is simply = (a0 ^ a1) || b1 || b0
 *
 * The differences above from [1] seem to improve the efficiency of evaulation
 * and they make the code easier to analyze.
 */

void ec_gen_table_2(const ecpt &a, const ecpt &b, ecpt TABLE[8]) {
	ecpt bn;
	ec_neg(b, bn);

	// P[4] = a
	ec_set(a, TABLE[4]);

	// P[5] = a - b
	ufe t2b;
	ec_add(a, bn, TABLE[5], true, true, true, t2b);

	// P[7] = a + b
	ec_add(a, b, TABLE[7], true, true, true, t2b);

	// P[6] = a + 2b
	ec_add(TABLE[7], b, TABLE[6], true, true, true, t2b);

	ecpt a2;
	ec_dbl(a, a2, true, t2b);

	// P[0] = 3a
	ec_add(a2, a, TABLE[0], true, false, true, t2b);

	// P[1] = 3a + b
	ec_add(TABLE[0], b, TABLE[1], true, true, true, t2b);

	// P[2] = 3a + 2b
	ec_add(TABLE[1], b, TABLE[2], true, true, true, t2b);

	// P[3] = 3a + 3b
	ec_add(TABLE[2], b, TABLE[3], true, true, true, t2b);
}

/*
 * Table index is simply = (a0 ^ a1) || b1 || b0
 */

static CAT_INLINE void ec_table_select_2(const ecpt *table, const ufp &a, const ufp &b, const int index, ecpt &r) {
	u32 bits = (u32)(a.w >> index);
	u32 k = ((bits ^ (bits >> 1)) & 1) << 2;
	k |= (u32)(b.w >> index) & 3;

	for (int ii = 0; ii < 8; ++ii) {
		// Generate a mask that is -1 if ii == index, else 0
		const u128 mask = ec_gen_mask(ii, k);

		// Add in the masked table entry
		ec_set_mask(table[ii], mask, r);
	}

	ec_cond_neg(((bits >> 1) & 1) ^ 1, r);
}

/*
 * Multiplication by variable base point
 *
 * Point must be in extended projective coordinates with T and Z values.
 *
 * The resulting point has undefined T and Z values, so must be expanded with
 * ec_expand() before using as input to other math functions.
 */

// R = kP
void ec_mul(const u64 k[4], const ecpt &P0, ecpt &R) {
	// Decompose scalar into subscalars
	ufp a, b;
	s32 asign, bsign;
	gls_decompose(k, asign, a, bsign, b);

	// Q = endomorphism of P
	ecpt Q;
	gls_morph(P0.x, P0.y, Q.x, Q.y);
	ec_expand(Q);

	// Set base point signs
	ecpt P;
	ec_set(P0, P);
	ec_cond_neg(asign, P);
	ec_cond_neg(bsign, Q);

	// Precompute multiplication table
	ecpt table[8];
	ec_gen_table_2(P, Q, table);

	// Recode subscalars
	u32 recode_bit = ec_recode_scalars_2(a, b, 128);

	// Initialize working point
	ecpt X;
	ec_table_select_2(table, a, b, 126, X);

	ufe t2b;
	for (int ii = 124; ii >= 0; ii -= 2) {
		ecpt T;
		ec_table_select_2(table, a, b, ii, T);

		ec_dbl(X, X, false, t2b);
		ec_dbl(X, X, false, t2b);
		ec_add(X, T, X, false, false, false, t2b);
	}

	// If bit == 1, X <- X + P (inverted logic from [1])
	ec_cond_add(recode_bit, X, P, X, true, false, t2b);

	// Compute affine coordinates in R
	ec_affine(X, false, R);
}


//// Constant-time Simultaneous Multiplication

/*
 * Precomputed table generation
 *
 * Using GLV-SAC Precomputation with m=4 [1], assuming window size of 1 bit
 */

void ec_gen_table_4(const ecpt &a, const ecpt &b, const ecpt &c, const ecpt &d, ecpt TABLE[8]) {
	// P[0] = a
	ec_set(a, TABLE[0]);

	// P[1] = a + b
	ufe t2b;
	ec_add(a, b, TABLE[1], true, true, true, t2b);

	// P[2] = a + c
	ec_add(a, c, TABLE[2], true, true, true, t2b);

	// P[3] = a + b + c
	ec_add(TABLE[1], c, TABLE[3], true, true, true, t2b);

	// P[4] = a + d
	ec_add(a, d, TABLE[4], true, true, true, t2b);

	// P[5] = a + b + d
	ec_add(TABLE[1], d, TABLE[5], true, true, true, t2b);

	// P[6] = a + c + d
	ec_add(TABLE[2], d, TABLE[6], true, true, true, t2b);

	// P[7] = a + b + c + d
	ec_add(TABLE[3], d, TABLE[7], true, true, true, t2b);
}

/*
 * Constant-time table selection for m=4
 */

static CAT_INLINE void ec_table_select_4(const ecpt *table, const ufp &a, const ufp &b, const ufp &c, const ufp &d, const int index, ecpt &r) {
	int k = ((u32)(b.w >> index) & 1);
	k |= ((u32)(c.w >> index) & 1) << 1;
	k |= ((u32)(d.w >> index) & 1) << 2;

	const int TABLE_SIZE = 8;
	for (int ii = 0; ii < TABLE_SIZE; ++ii) {
		// Generate a mask that is -1 if ii == index, else 0
		const u128 mask = ec_gen_mask(ii, k);

		// Add in the masked table entry
		ec_set_mask(table[ii], mask, r);
	}

	ec_cond_neg(((a.w >> index) & 1) ^ 1, r);
}

/*
 * Simultaneous multiplication by two variable base points
 *
 * Points must be in extended projective coordinates with T and Z values.
 *
 * The resulting point has undefined T and Z values, so must be expanded with
 * ec_expand() before using as input to other math functions.
 */

// R = aP + bQ
void ec_simul(const u64 a[4], const ecpt &P, const u64 b[4], const ecpt &Q, ecpt &R) {
	// Decompose scalar into subscalars
	ufp a0, a1, b0, b1;
	s32 a0sign, a1sign, b0sign, b1sign;
	gls_decompose(a, a0sign, a0, a1sign, a1);
	gls_decompose(b, b0sign, b0, b1sign, b1);

	// P1, Q1 = endomorphism points
	ecpt P1, Q1;
	gls_morph(P.x, P.y, P1.x, P1.y);
	ec_expand(P1);
	gls_morph(Q.x, Q.y, Q1.x, Q1.y);
	ec_expand(Q1);

	// Set base point signs
	ecpt P0, Q0;
	ec_set(P, P0);
	ec_set(Q, Q0);
	ec_cond_neg(a0sign, P0);
	ec_cond_neg(b0sign, Q0);
	ec_cond_neg(a1sign, P1);
	ec_cond_neg(b1sign, Q1);

	// Precompute multiplication table
	ecpt table[8];
	ec_gen_table_4(P0, P1, Q0, Q1, table);

	// Recode scalar
	u32 recode_bit = ec_recode_scalars_4(a0, a1, b0, b1, 127);

	// Initialize working point
	ecpt X;
	ec_table_select_4(table, a0, a1, b0, b1, 126, X);

	ufe t2b;
	for (int ii = 125; ii >= 0; --ii) {
		ecpt T;
		ec_table_select_4(table, a0, a1, b0, b1, ii, T);

		ec_dbl(X, X, false, t2b);
		ec_add(X, T, X, false, false, false, t2b);
	}

	// If bit == 1, X <- X + P (inverted logic from [1])
	ec_cond_add(recode_bit, X, P0, X, true, false, t2b);

	// Compute affine coordinates in R
	ec_affine(X, false, R);
}

