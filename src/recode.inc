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

#ifdef CAT_HAS_VECTOR_EXTENSIONS

		const vec_ecpt *tp = (const vec_ecpt *)table;
		vec_ecpt *rp = (vec_ecpt *)&r;

		for (int ii = 0; ii < 8; ++ii) {
			// Generate a mask that is -1 if ii == index, else 0
			const u64 mask = ec_gen_mask(ii, k);

			*rp ^= tp[ii] & mask;
		}

#else

		for (int ii = 0; ii < 8; ++ii) {
			// Generate a mask that is -1 if ii == index, else 0
			const u64 mask = ec_gen_mask(ii, k);

			ec_xor_mask(table[ii], mask, r);
		}

#endif
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

#ifdef CAT_HAS_VECTOR_EXTENSIONS

	const vec_ecpt *tp = (const vec_ecpt *)table;
	vec_ecpt *rp = (vec_ecpt *)&r;

	for (int ii = 0; ii < 8; ++ii) {
		// Generate a mask that is -1 if ii == index, else 0
		const u64 mask = ec_gen_mask(ii, k);

		*rp ^= tp[ii] & mask;
	}
	
#else

	for (int ii = 0; ii < 8; ++ii) {
		// Generate a mask that is -1 if ii == index, else 0
		const u64 mask = ec_gen_mask(ii, k);

		// Add in the masked table entry
		ec_xor_mask(table[ii], mask, r);
	}

#endif

	ec_cond_neg(((a.w >> index) & 1) ^ 1, r);
}

/*
 * Recode scalar to signed binary representation from [17]
 *
 * A good algorithm was hinted at in the paper, which can be summarized as
 * exploiting the group order to produce a signed encoding in the even case
 * by representing x - p, or x - 2p in the odd case.  So the resulting
 * encoding is the same length as the original bitstring and it does not
 * require any corrective adds/negates after evaluation.
 *
 * This turns out to be simpler than the comb approach from [1]
 * both for recoding and evaluation.
 *
 * Preconditions:
 *	k < q
 */

static CAT_INLINE void ec_recode_comb(const u64 k[4], u64 b[4]) {
	// Generate 64-bit mask from low bit
	const u64 lsb = k[0] & 1;
	const s64 mask = -(s64)lsb;

	// a = (k0 ? Q2 : Q1), where Q1 has an additional high bit set
	const u64 a3 = (lsb ^ 1) << 59;
	const u64 a0 = 0x18B24B8E27B0FF2DULL ^ (0x29D6DC9268D10176ULL & mask);
	const u64 a1 = 0x2CECF5F59F91BC16ULL ^ (0x75351E1EA0B2C43AULL & mask);

	// b = (k >> 1) + a
	u128 sum = (u128)((k[0] >> 1) | (k[1] << 63)) + a0;
	b[0] = (u64)sum;
	sum = ((u128)((k[1] >> 1) | (k[2] << 63)) + a1) + (u64)(sum >> 64);
	b[1] = (u64)sum;
	sum = (u128)((k[2] >> 1) | (k[3] << 63)) + (u64)(sum >> 64);
	b[2] = (u64)sum;
	sum = ((u128)(k[3] >> 1) + a3) + (u64)(sum >> 64);
	b[3] = (u64)sum;
}

/*
 * LSB-Set Comb Method Scalar Recoding [1] for w=7, v=6
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
 * Number of ECADDs and ECDBLs = d - 1, e - 1, respectively.
 *
 * Each step during evaulation needs to access all of the table entries,
 * which takes time.  So the total number of entries needs to be tuned along
 * with w,v.
 *
 * Memory accesses cost 1530 cycles for one 128-entry table lookup.
 * ec_dbl cost 450 cycles.
 * ec_add cost 570 cycles.
 *
 * Rough cost = (e - 1) * 570 + (d - 1) * 450 + (1530/128)*2^(w-1)*v * d
 */

/*
	Modern ARM and Intel processors all have L1 data cache size of 32 KB,
	so going much above 16 KB seems like a bad idea.

	v	w	e	d	cost	table size	table bytes

	These are too big for L1 data cache of 32 KB:

	6	8	6	36	77280	768	49152
	5	8	7	35	75630	640	40960
	2	9	14	28	106920	512	32768
	4	8	8	32	69780	512	32768
	8	7	5	40	54630	512	32768
	7	7	6	42	57750	448	28672
	3	8	11	33	73230	384	24576
	6	7	6	36	49740	384	24576
	5	7	8	40	55980	320	20480

	Interesting region:

	1	9	28	28	113220	256	16384
	4	7	9	36	51090	256	16384 <- ec_mul_gen: 59805 cycles 22 usec
	8	6	6	48	47400	256	16384 <- Longa's parameters [1] on this curve
	7	6	6	42	41685	224	14336 <- ec_mul_gen: 51402 cycles 19 usec
	3	7	12	36	52440	192	12288 <- ec_mul_gen: 60048 cycles 22 usec
	6	6	7	42	42135	192	12288 <- ec_mul_gen: 56979 cycles 21 usec
	5	6	9	45	45892.5	160	10240
	2	7	18	36	55140	128	8192 <- ec_mul_gen: 58890 cycles 22 usec
	4	6	11	44	45840	128	8192 <- Longa's parameters [1] on this curve
	8	5	7	56	44760	128	8192
	7	5	8	56	45210	112	7168
	3	6	14	42	45285	96	6144 <- ec_mul_gen: 53517 cycles 19 usec
	6	5	9	54	44137.5	96	6144
	5	5	11	55	45799	80	5120 <- Hamburg's parameters [17] on this curve
	1	7	36	36	63240	64	4096
	2	6	21	42	48435	64	4096 <- ec_mul_gen: 56442 cycles 20 usec
	4	5	13	52	44415	64	4096 <- Longa's parameters [1] on this curve
	3	5	17	51	45454	48	3072
	1	6	42	42	57885	32	2048
	2	5	26	52	50265	32	2048
	1	5	51	51	60754	16	1024

	I tested a number of configurations, neglecting v=5 or w=5 since it does
	not evenly divide t=252.  The cost function was not entirely accurate.
 */

// Selected parameters for ec_mul_gen:
static const int MG_l = 252;
static const int MG_t = 6; // number of comb teeth
static const int MG_n = 7; // number of combs
static const int MG_s = MG_l / (MG_n * MG_t); // = ceil(t/wv), comb teeth spacing
static const int MG_w = MG_s * MG_t; // width of a comb
static const int MG_e = 1 << (MG_t - 1); // number of comb table entries

static void ec_table_select_comb_gen(const u64 b[4], const int ii, ecpt r[MG_n]) {
	// For each comb,
	for (int comb = 0, comb_bit = ii; comb < MG_n; ++comb, comb_bit += MG_w) {
		u32 bits = 0, offset = comb_bit;
		for (int tooth = 0; tooth < (MG_t-1); ++tooth) {
			bits |= ((u32)(b[offset >> 6] >> (offset & 63)) & 1) << tooth;
			offset += MG_s;
		}
		u32 sign = ((u32)(b[offset >> 6] >> (offset & 63)) & 1) ^ 1;
		bits = (bits ^ -(s32)sign) & (MG_e - 1);

		ecpt &p = r[comb];

		ec_zero(p);

		const ecpt_affine *table = GEN_TABLE[comb];

#ifdef CAT_HAS_VECTOR_EXTENSIONS

		const vec_ecpt_affine *tp = (const vec_ecpt_affine *)table;
		vec_ecpt_affine *rp = (vec_ecpt_affine *)&p;

		for (int jj = 0; jj < MG_e; ++jj) {
			// Generate a mask that is -1 if jj == index, else 0
			const u64 mask = ec_gen_mask(jj, bits);

			*rp ^= tp[jj] & mask;
		}
	
#else

		for (int jj = 0; jj < MG_e; ++jj) {
			// Generate a mask that is -1 if jj == index, else 0
			const u64 mask = ec_gen_mask(jj, bits);

			// Add in the masked table entry
			ec_xor_mask_affine(table[jj], mask, p);
		}

#endif

		// Reconstruct T
		fe_mul(p.x, p.y, p.t);

		// Apply sign bit
		ec_cond_neg(sign ^ 1, p);
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

static CAT_INLINE u32 comb_bit_81(const u64 b[4], const int wp, const int ep) {
	// K(w', v', e') = b_(d * w' + e * v' + e'), v' = 0
	const u32 jj = (wp << 5) + ep;

	return (u32)(b[jj >> 6] >> (jj & 63)) & 1;
}

// NOTE: Not constant time because it does not need to be for ec_simul_gen
static void ec_table_select_comb_81(const u64 b[4], const int ii, ecpt &p) {
	// D(v', e') = K(w-1, v', e') || K(w-2, v', e') || ... || K(1, v', e')
	// s(v', e') = K(0, v', e')

	u32 d = comb_bit_81(b, 7, ii) << 6;
	d |= comb_bit_81(b, 6, ii) << 5;
	d |= comb_bit_81(b, 5, ii) << 4;
	d |= comb_bit_81(b, 4, ii) << 3;
	d |= comb_bit_81(b, 3, ii) << 2;
	d |= comb_bit_81(b, 2, ii) << 1;
	d |= comb_bit_81(b, 1, ii);
	const u32 s = comb_bit_81(b, 0, ii) ^ 1;

	p.x = SIMUL_GEN_TABLE[d].x;
	p.y = SIMUL_GEN_TABLE[d].y;
	p.t = SIMUL_GEN_TABLE[d].t;

	// Flip recode_lsb sign here rather than at the end to interleave easier
	if (s) {
		ec_neg(p, p);
	}
}

