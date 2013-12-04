// Elliptic curve point multiplication scalar recoding and table generation

#include "precomp.inc"

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

static CAT_INLINE void ec_gen_table_2(const ecpt &a, const ecpt &b, ecpt TABLE[8]) {
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
 * GLV-SAC Scalar Recoding Algorithm for m=2 [1]
 *
 * Returns low bit of 'a'
 */

static CAT_INLINE u32 ec_recode_scalars_2(ufp &a, ufp &b, const int len) {
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
 * Table index is simply = (a0 ^ a1) || b1 || b0
 */

static void ec_table_select_2(const ecpt *table, const ufp &a, const ufp &b, const int index, const bool constant_time, ecpt &r) {
	u32 bits = (u32)(a.w >> index);
	u32 k = ((bits ^ (bits >> 1)) & 1) << 2;
	k |= (u32)(b.w >> index) & 3;

	// If constant time requested,
	if (constant_time) {
		ec_zero(r);

		for (int ii = 0; ii < 8; ++ii) {
			// Generate a mask that is -1 if ii == index, else 0
			const u64 mask = ec_gen_mask(ii, k);

			// Add in the masked table entry
			ec_xor_mask(table[ii], mask, r);
		}
	} else {
		ec_set(table[k], r);
	}

	ec_cond_neg(((bits >> 1) & 1) ^ 1, r);
}

/*
 * Precomputed table generation
 *
 * Using GLV-SAC Precomputation with m=4 [1], assuming window size of 1 bit
 */

static void ec_gen_table_4(const ecpt &a, const ecpt &b, const ecpt &c, const ecpt &d, ecpt TABLE[8]) {
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
 * GLV-SAC Scalar Recoding Algorithm for m=4 [1]
 *
 * Returns low bit of 'a'
 */

static CAT_INLINE u32 ec_recode_scalars_4(ufp &a, ufp &b, ufp &c, ufp &d, const int len) {
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

/*
 * Constant-time table selection for m=4
 */

static void ec_table_select_4(const ecpt *table, const ufp &a, const ufp &b, const ufp &c, const ufp &d, const int index, ecpt &r) {
	int k = ((u32)(b.w >> index) & 1);
	k |= ((u32)(c.w >> index) & 1) << 1;
	k |= ((u32)(d.w >> index) & 1) << 2;

	ec_zero(r);

	const int TABLE_SIZE = 8;
	for (int ii = 0; ii < TABLE_SIZE; ++ii) {
		// Generate a mask that is -1 if ii == index, else 0
		const u64 mask = ec_gen_mask(ii, k);

		// Add in the masked table entry
		ec_xor_mask(table[ii], mask, r);
	}

	ec_cond_neg(((a.w >> index) & 1) ^ 1, r);
}

/*
 * LSB-Set Comb Method Scalar Recoding [1] for w=7, v=2
 *
 * The algorithm is tuned with ECADD = 1.64 * ECDBL in cycles.
 *
 * t = 252 bits for input scalars
 * w = window size in bits
 * v = number of tables
 * e = roundup(t / wv)
 * d = e * v
 * l = d * w
 *
 * The parameters w,v are tunable.  The number of table entries:
 *	= v * 2^(w - 1)
 *
 * Number of ECADDs and ECDBLs = e - 1, d - 1, respectively.
 *
 * Constant-time table lookups are expensive and proportional to
 * the number of overall entries.  After some experimentation,
 * 256 is too many table entries.  Halving it to 128 entries is
 * a good trade-off.
 *
 * Optimizing for operation counts, choosing values of v,w that
 * yield tables of 128 entries:
 *
 * v,w -> e,d   -> effective ECDBLs
 * 1,8 -> 32,32 -> 81.84
 * 2,7 -> 18,36 -> 74.4 <- best option
 * 4,6 -> 11,44 -> 80.52
 */

static CAT_INLINE u32 ec_recode_scalar_comb(const u64 k[4], u64 b[4]) {
	const int t = 252;
	const int w = 7;
	const int v = 2;
	const int e = t / (w * v); // t / wv
	const int d = e * v; // ev
	const int l = d * w; // dw

	// If k0 == 0, b = q - k (and return 1), else b = k (and return 0)

	const u32 lsb = (u32)k[0] & 1;
	const u64 mask = (s64)0 - lsb;

	neg_mod_q(k, b);

	b[0] ^= (k[0] ^ b[0]) & mask;
	b[1] ^= (k[1] ^ b[1]) & mask;
	b[2] ^= (k[2] ^ b[2]) & mask;
	b[3] ^= (k[3] ^ b[3]) & mask;

	// Recode scalar:

	const u64 d_bit = (u64)1 << (d - 1);
	const u64 low_mask = d_bit - 1;

	// for bits 0..(d-1), 0 => -1, 1 => +1
	b[0] = (b[0] | (low_mask | d_bit)) ^ (d_bit | ((b[0] >> 1) & low_mask));

	for (int i = d; i < l; ++i) {
		u32 b_imd = (u32)(b[0] >> (i % d));
		u32 b_i = (u32)(b[i >> 6] >> (i & 63));
		u32 bit = b_imd & b_i & 1;

		const int j = i + 1;
		u64 t[4] = {0};
		t[j >> 6] |= (u64)bit << (j & 63);

		u128 sum = (u128)b[0] + t[0];
		b[0] = (u64)sum;
		sum = ((u128)b[1] + t[1]) + (u64)(sum >> 64);
		b[1] = (u64)sum;
		sum = ((u128)b[2] + t[2]) + (u64)(sum >> 64);
		b[2] = (u64)sum;
		sum = ((u128)b[3] + t[3]) + (u64)(sum >> 64);
		b[3] = (u64)sum;
	}

	return lsb ^ 1;
}

static CAT_INLINE u32 comb_bit(const u64 b[4], const int wp, const int vp, const int ep) {
	// K(w', v', e') = b_(d * w' + e * v' + e')
	u32 jj = (wp * 36) + (vp * 18) + ep;

	return (u32)(b[jj >> 6] >> (jj & 63)) & 1;
}

static void ec_table_select_comb(const u64 b[4], const int ii, const bool constant_time, ecpt &p1, ecpt &p2) {
	// DCK(v', e') = K(w-1, v', e') || K(w-2, v', e') || ... || K(1, v', e')
	// s(v', e') = K(0, v', e')

	// Select table entry 
	// p1 = s(0, ii) * tables[DCK(0, ii)][0]
	// p2 = s(1, ii) * tables[DCK(1, ii)][1]

	u32 d_0;
	d_0 = comb_bit(b, 6, 0, ii) << 5;
	d_0 |= comb_bit(b, 5, 0, ii) << 4;
	d_0 |= comb_bit(b, 4, 0, ii) << 3;
	d_0 |= comb_bit(b, 3, 0, ii) << 2;
	d_0 |= comb_bit(b, 2, 0, ii) << 1;
	d_0 |= comb_bit(b, 1, 0, ii);
	u32 s_0 = comb_bit(b, 0, 0, ii);

	if (constant_time) {
		ec_zero(p1);
		for (int jj = 0; jj < 64; ++jj) {
			// Generate a mask that is -1 if jj == index, else 0
			const u64 mask = ec_gen_mask(jj, d_0);

			// Add in the masked table entry
			ec_xor_mask_affine(GEN_TABLE_0[jj], mask, p1);
		}
		fe_mul(p1.x, p1.y, p1.t);
		ec_cond_neg(s_0, p1);
	} else {
		ec_expand(GEN_TABLE_0[d_0], p1);
		if (s_0) {
			ec_neg(p1, p1);
		}
	}

	u32 d_1;
	d_1 = comb_bit(b, 6, 1, ii) << 5;
	d_1 |= comb_bit(b, 5, 1, ii) << 4;
	d_1 |= comb_bit(b, 4, 1, ii) << 3;
	d_1 |= comb_bit(b, 3, 1, ii) << 2;
	d_1 |= comb_bit(b, 2, 1, ii) << 1;
	d_1 |= comb_bit(b, 1, 1, ii);
	u32 s_1 = comb_bit(b, 0, 1, ii);

	if (constant_time) {
		ec_zero(p2);
		for (int jj = 0; jj < 64; ++jj) {
			// Generate a mask that is -1 if jj == index, else 0
			const u64 mask = ec_gen_mask(jj, d_1);

			// Add in the masked table entry
			ec_xor_mask_affine(GEN_TABLE_1[jj], mask, p2);
		}
		fe_mul(p2.x, p2.y, p2.t);
		ec_cond_neg(s_1, p2);
	} else {
		ec_expand(GEN_TABLE_1[d_1], p2);
		if (s_1) {
			ec_neg(p2, p2);
		}
	}
}

/*
 * LSB-set Scalar Recoding [1] with w=8, v=1
 *
 * This function is useful for EdDSA signature validation, so it is
 * interesting to optimize for this case.
 *
 * Interleaving the ECADDs for ec_mul with those from ec_mul_gen is
 * a straight-forward approach.  We want the ec_mul_gen table to
 * stay at 128 points since that is the optimal memory access time
 * trade-off.  But, there is no need to use multiple tables since
 * the ECDBLs need to be performed *anyway* for the ec_mul ops,
 * so the ECDBLs are sort-of "free."  So the optimal choice for
 * table construction is a little different from the ec_mul_gen
 * case and we need a new table for w = 8, v = 1.  Since 8 does
 * not evenly divide 252, it is not necessary to do the final
 * correction step addition which simplifies the algorithm a bit.
 *
 * For this tuning, ec_mul_gen ECADDs = 32.
 *
 * Since ECDBL/ECADD ops are linear, it is possible to interleave
 * ec_mul_gen and ec_mul even though the number of ECDBL for each
 * is different.  Introducing ECADDs for ec_mul_gen near the end
 * of the evaluation loop of ec_mul still exhibits a regular
 * pattern and will just require another 32 ECADDs.  The final
 * conditional negation from ec_mul_gen can be merged into the
 * ECADDs by inverting the sign of each added point instead to
 * avoid messing with the interleaving.
 *
 * So overall the cost should be about the same as one ec_mul
 * with just 32 extra ECADDs from table lookups, which falls
 * about mid-way between ec_mul and ec_simul for performance.
 */

static u32 ec_recode_scalar_comb1(const u64 k[4], u64 b[4]) {
	//const int t = 252;
	const int w = 8;
	const int v = 1;
	const int e = 32; // t / wv
	const int d = e * v;
	const int l = d * w;

	// If k0 == 0, b = q - k (and return 1), else b = k (and return 0)

	const u32 lsb = (u32)k[0] & 1;
	const u64 mask = (s64)0 - lsb;

	neg_mod_q(k, b);

	b[0] ^= (k[0] ^ b[0]) & mask;
	b[1] ^= (k[1] ^ b[1]) & mask;
	b[2] ^= (k[2] ^ b[2]) & mask;
	b[3] ^= (k[3] ^ b[3]) & mask;

	// Recode scalar:

	const u64 d_bit = (u64)1 << (d - 1);
	const u64 low_mask = d_bit - 1;

	// for bits 0..(d-1), 0 => -1, 1 => +1
	b[0] = (b[0] | (low_mask | d_bit)) ^ (d_bit | ((b[0] >> 1) & low_mask));

	for (int i = d; i < l; ++i) {
		u32 b_imd = (u32)(b[0] >> (i & 31));
		u32 b_i = (u32)(b[i >> 6] >> (i & 63));
		u32 bit = b_imd & b_i & 1;

		const int j = i + 1;
		u64 t[4] = {0};
		t[j >> 6] |= (u64)bit << (j & 63);

		u128 sum = (u128)b[0] + t[0];
		b[0] = (u64)sum;
		sum = ((u128)b[1] + t[1]) + (u64)(sum >> 64);
		b[1] = (u64)sum;
		sum = ((u128)b[2] + t[2]) + (u64)(sum >> 64);
		b[2] = (u64)sum;
		sum = ((u128)b[3] + t[3]) + (u64)(sum >> 64);
		b[3] = (u64)sum;
	}

	return lsb ^ 1;
}

static CAT_INLINE u32 comb_bit1(const u64 b[4], const int wp, const int ep) {
	// K(w', v', e') = b_(d * w' + e * v' + e')
	const u32 jj = (wp << 5) + ep;
	return (u32)(b[jj >> 6] >> (jj & 63)) & 1;
}

// NOTE: Not constant time because it does not need to be for ec_simul_gen
static void ec_table_select_comb1(const u32 recode_lsb, const u64 b[4], const int ii, ecpt &p1) {
	// DCK(v', e') = K(w-1, v', e') || K(w-2, v', e') || ... || K(1, v', e')
	// s(v', e') = K(0, v', e')

	// Select table entry 
	// p1 = s(0, ii) * tables[DCK(0, ii)][0]
	// p2 = s(1, ii) * tables[DCK(1, ii)][1]

	u32 d_0;
	d_0 = comb_bit1(b, 7, ii) << 6;
	d_0 |= comb_bit1(b, 6, ii) << 5;
	d_0 |= comb_bit1(b, 5, ii) << 4;
	d_0 |= comb_bit1(b, 4, ii) << 3;
	d_0 |= comb_bit1(b, 3, ii) << 2;
	d_0 |= comb_bit1(b, 2, ii) << 1;
	d_0 |= comb_bit1(b, 1, ii);
	u32 s_0 = comb_bit1(b, 0, ii);

	ec_expand(SIMUL_GEN_TABLE[d_0], p1);

	// Flip sign here rather than at the end to interleave more easily
	if (s_0 ^ recode_lsb) {
		ec_neg(p1, p1);
	}
}
