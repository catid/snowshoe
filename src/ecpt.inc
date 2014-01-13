// Endomorphism
#include "endo.inc"
#include "ecpt.hpp"

/*
 * "Generator point": The smallest X-coordinate point of q torsion
 *
 * Gx: 0 * i + 15
 * Gy: 0x6E848B46758BA443DD9869FE923191B0 * i + 0x7869C919DD649B4C36D073DADE2014AB
 * Gt = Gx * Gy
 * = 0x79C42920E32E9FF9FBEE35EA90E7895C * i + 0xE32C883F8E519773636C9D303E13613
 * Gz = 1
 */

static const ufe EC_GX = {
	{{15}}, {{0}}
};

static const ufe EC_GY = {
	{{	// a
		0x36D073DADE2014ABULL,
		0x7869C919DD649B4CULL
	}},
	{{	// b
		0xDD9869FE923191B0ULL,
		0x6E848B46758BA443ULL
	}}
};

static const ufe EC_GT = {
	{{	// a
		0x3636C9D303E13613ULL,
		0x0E32C883F8E51977ULL
	}},
	{{	// b
		0xFBEE35EA90E7895CULL,
		0x79C42920E32E9FF9ULL
	}}
};

static const ufe EC_1 = {
	{{ 1, 0 }},
	{{ 0, 0 }}
};


static const ecpt EC_G = {
	EC_GX, EC_GY, EC_GT, EC_1
};

/*
 * Endomorphism of generator point
 *
 * Calculated as (wx*conj(x), conj(y)) = (x2, y2)
 *
 * x2 = 13933292142614162712970301776591149603*i
 *      + 92037237801541697222328802746237627665
 * y2 = 23237895901517042326583462276207898191*i
 * 		+ 160056629477888271617002104343968617643
 *
 * x2 = 0xA7B74573CBFDEDC58A1315F74055A23 * i
 *	  + 0x453DBA2B9E5FEF6E2C5098AFBA02AD11
 * y2 = 0x117B74B98A745BBC226796016DCE6E4F * i
 *	  + 0x7869C919DD649B4C36D073DADE2014AB
 *
 * t2 = x2 * y2
 * = 0x7DDDF1B7EDE526FD7E31219AE8A9CF59 * i + 0x6CF2049A858BC394ED7E4F54BAAE9EA3
 */

static const ufe EC_EGX = {
	{{	// a
		0x2C5098AFBA02AD11ULL,
		0x453DBA2B9E5FEF6EULL
	}},
	{{	// b
		0x58A1315F74055A23ULL,
		0x0A7B74573CBFDEDCULL
	}}
};

static const ufe EC_EGY = {
	{{	// a
		0x36D073DADE2014ABULL,
		0x7869C919DD649B4CULL
	}},
	{{	// b
		0x226796016DCE6E4FULL,
		0x117B74B98A745BBCULL
	}}
};

static const ufe EC_EGT = {
	{{	// a
		0xED7E4F54BAAE9EA3ULL,
		0x6CF2049A858BC394ULL
	}},
	{{	// b
		0x7E31219AE8A9CF59ULL,
		0x7DDDF1B7EDE526FDULL
	}}
};

static const ecpt EC_EG = {
	EC_EGX, EC_EGY, EC_EGT, EC_1
};

// Load (x,y) from endian-neutral data bytes (64)
static void ec_load_xy(const u8 *a, ecpt_affine &r) {
	fe_load(a, r.x);
	fe_load(a + 32, r.y);
}

// Save (x,y) to endian-neutral data bytes (64)
static void ec_save_xy(const ecpt_affine &a, u8 *r) {
	fe_save(a.x, r);
	fe_save(a.y, r + 32);
}

// r = 0
static CAT_INLINE void ec_zero(ecpt &r) {
	fe_zero(r.x);
	fe_zero(r.y);
	fe_zero(r.t);
	fe_zero(r.z);
}

// r = identity element s.t. a + r = a
static CAT_INLINE void ec_identity(ecpt &r) {
	fe_zero(r.x);
	fe_set_smallk(1, r.y);
	fe_zero(r.t);
	fe_set_smallk(1, r.z);
}

// r = -a
static CAT_INLINE void ec_neg(const ecpt &a, ecpt &r) {
	// -(X : Y : T : Z) = (-X : Y : -T : Z)
	fe_neg(a.x, r.x);
	fe_set(a.y, r.y);
	fe_neg(a.t, r.t);
	fe_set(a.z, r.z);
}

// r = -a
static CAT_INLINE void ec_neg_affine(const ecpt_affine &a, ecpt_affine &r) {
	// -(X, Y) = (-X, Y)
	fe_neg(a.x, r.x);
	fe_set(a.y, r.y);
}

// r = (mask==-1 ? r : -r)
static CAT_INLINE void ec_neg_mask_inplace(const u64 mask, ecpt &r) {
	fe_neg_mask_inplace(mask, r.x);
	fe_neg_mask_inplace(mask, r.t);
}

// r = (mask==-1 ? p : -p)
static CAT_INLINE void ec_neg_mask(const u64 mask, const ecpt &p, ecpt &r) {
	fe_neg_mask(mask, p.x, r.x);
	fe_set(p.y, r.y);
	fe_neg_mask(mask, p.t, r.t);
	fe_set(p.z, r.z);
}

// r = (mask==-1 ? r : -r)
static CAT_INLINE void ec_neg_mask_affine(const u64 mask, ecpt_affine &r) {
	fe_neg_mask(mask, r.x, r.x);
}

// r = a
static CAT_INLINE void ec_set(const ecpt &a, ecpt &r) {
	fe_set(a.x, r.x);
	fe_set(a.y, r.y);
	fe_set(a.t, r.t);
	fe_set(a.z, r.z);
}

// Expand from affine coordinates to extended
static CAT_INLINE void ec_expand(const ecpt_affine &a, ecpt &r) {
	fe_set(a.x, r.x);
	fe_set(a.y, r.y);

	// t = xy, z = 1
	fe_mul(r.x, r.y, r.t);
	fe_set_smallk(1, r.z);
}

// r = (mask == -1) ? a : r
static CAT_INLINE void ec_set_mask(const ecpt &a, const u64 mask, ecpt &r) {
	fe_set_mask(a.x, mask, r.x);
	fe_set_mask(a.y, mask, r.y);
	fe_set_mask(a.t, mask, r.t);
	fe_set_mask(a.z, mask, r.z);
}

// r ^= a & mask
static CAT_INLINE void ec_xor_mask(const ecpt &a, const u64 mask, ecpt &r) {

#ifdef CAT_HAS_VECTOR_EXTENSIONS

	const vec_ecpt *tp = (const vec_ecpt *)&a;
	vec_ecpt *rp = (vec_ecpt *)&r;

	*rp ^= *tp & mask;
	
#else

	fe_xor_mask(a.x, mask, r.x);
	fe_xor_mask(a.y, mask, r.y);
	fe_xor_mask(a.t, mask, r.t);
	fe_xor_mask(a.z, mask, r.z);

#endif

}

// r ^= a & mask
static CAT_INLINE void ec_xor_mask_affine(const ecpt_affine &a, const u64 mask, ecpt &r) {

#ifdef CAT_HAS_VECTOR_EXTENSIONS

	const vec_ecpt_affine *tp = (const vec_ecpt_affine *)&a;
	vec_ecpt_affine *rp = (vec_ecpt_affine *)&r;

	*rp ^= *tp & mask;
	
#else

	fe_xor_mask(a.x, mask, r.x);
	fe_xor_mask(a.y, mask, r.y);

#endif

}

/*
 * Extended Twisted Edwards Doubling
 *
 * (x2, y2, t2, z2) = 2 * (x, y, z), where t2 = x2 * y2 / z2
 *
 * This doubling formula is based on the dedicated doubling formula from [5].
 * Optimized for Instruction-Level Paralellism (ILP) and low register usage.
 *
 * The original formula appears to produce all negative results, which is okay.
 * It was adjusted to produce the same as the math expressions in the paper,
 * which works out better for my 'u'.
 *
 * There are alternative formulae in [3] that do not seem preferable for my Fp.
 * The EFD website also has a z1=1 entry that is more efficient, but would
 * require implementing a whole new ec_dbl function so does not seem worth it.
 *
 * This doubling formula produces a split value for T similar to Hamburg [17].
 * The value of T can be reconstructed after ec_dbl by multiplying r.t and t2b.
 * This lazier calculation is preferable to generating T in this function so
 * that a regular execution pattern can be used to defeat timing attack.
 */

// r = 2p
static void ec_dbl(const ecpt &p, ecpt &r, const bool z_one, ufe &t2b) {
	// Uses 4S 3M 7A 1U when calc_t=false, z_one=false
	// z_one=true: -1S -1A

	// t <- x + y
	fe_add(p.x, p.y, r.t);

	// z <- z^2
	if (!z_one) {
		fe_sqr(p.z, r.z);
	}

	// x <- x^2
	fe_sqr(p.x, r.x);

	// t <- t^2
	fe_sqr(r.t, r.t);

	// y <- y^2
	fe_sqr(p.y, r.y);

	// t <- t - x = (x + y)^2 - x^2
	fe_sub(r.t, r.x, r.t);

	// x <- u * x = u * x^2
	fe_mul_u(r.x, r.x);

	// z <- z + z = 2 * z^2
	if (!z_one) {
		fe_add(r.z, r.z, r.z);
	} else {
		fe_set_smallk(2, r.z);
	}

	// w <- y - x = y^2 - u * x^2
	ufe w;
	fe_sub(r.y, r.x, w);

	// t <- t - y = (x + y)^2 - x^2 - y^2 = 2 * x * y
	fe_sub(r.t, r.y, r.t);

	// z <- z - w = 2 * z^2 - y^2 + u * x^2
	fe_sub(r.z, w, r.z);

	// t2b <- y + x = y^2 + u * x^2
	fe_add(r.x, r.y, t2b);

	// x2 <- t * z = 2 * x * y * (2 * z^2 - y^2 + u * x^2)
	fe_mul(r.t, r.z, r.x);

	// t2 <- t * t2b = 2 * x * y * (y^2 + u * x^2)
	// Not performed here: Can be reconstructed from t2 = r.t * t2b later

	// y2 <- w * t2b = (y^2 - u * x^2) * (y^2 + u * x^2)
	fe_mul(w, t2b, r.y);

	// z2 <- w * z = (y^2 - u * x^2) * (2 * z^2 - y^2 + u * x^2)
	fe_mul(w, r.z, r.z);
}

/*
 * Extended Twisted Edwards Unified Point Addition
 *
 * There are a number of faster addition laws available for Twisted Edwards
 * elliptic curves [5].  However the robustness of a complete, unified law is
 * hard to argue against.  Especially since the evaluation stage of the point
 * multiplication code is only executing one ec_add for every 2 ec_dbl calls
 * for variable-base scalar multiplication (most important server-side math),
 * and faster addition laws only reduce the operation count by ~1 M at best.
 *
 * Optimized for Instruction-Level Paralellism (ILP) and low register usage.
 *
 * Both input points are in extended projective coordinates:
 * (X3, Y3, T3, Z3) = (X1, Y1, T1, Z1) + (X2, Y2, T2, Z2)
 *
 * This addition formula produces a split value for T similar to Hamburg [17].
 * The value of T can be reconstructed after ec_add by multiplying r.t and t2b.
 * This lazier calculation is preferable to generating T in this function so
 * that a regular execution pattern can be used to defeat timing attack.
 *
 * If in_precomp_t1 is true:
 *	The t2b parameter will be ignored.  The p1.t field is assumed to be T1.
 * Else:
 *	The t2b parameter will be multiplied by p1.t to recover T1.
 *
 * If out_precomp_t3 is true:
 *	The t2b parameter is undefined.  The r.t field contains T1.
 * Else:
 *	This function generates r.t and t2b as partial products.
 */

// r = p1 + p2
static void ec_add(const ecpt &p1, const ecpt &p2, ecpt &r, const bool z2_one, const bool in_precomp_t1, const bool out_precomp_t3, ufe &t2b) {
	// Uses: 9M 7A 1D 2U with all flags false
	// z2_one=true: -1M
	// in_precomp_t1=true: -1M
	// out_precomp_t1=true: +1M

	// If t1 is not precomputed in p1.t,
	if (!in_precomp_t1) {
		// t1 <- t1 * t2b
		fe_mul(p1.t, t2b, r.t);
	}

	// t2b <- (x1 + y1)
	fe_add(p1.x, p1.y, t2b);

	// w1 <- (x2 + y2)
	ufe w1;
	fe_add(p2.x, p2.y, w1);

	// w2 <- t1 * t2
	ufe w2;
	if (in_precomp_t1) {
		fe_mul(p1.t, p2.t, w2);
	} else {
		fe_mul(r.t, p2.t, w2);
	}

	// t2b <- t2b * w1 = (x1 + y1) (x2 + y2)
	fe_mul(t2b, w1, t2b);

	// t <- x1 * x2
	fe_mul(p1.x, p2.x, r.t);

	// y <- y1 * y2
	fe_mul(p1.y, p2.y, r.y);

	// w2 <- u * w2 = u * t1 * t2
	fe_mul_u(w2, w2);

	// t2b <- t2b - t = (x1 + y1) (x2 + y2) - x1 * x2
	fe_sub(t2b, r.t, t2b);

	// w2 <- d * w2 = u * d * t1 * t2
	fe_mul_smallk(w2, EC_D, w2);

	// t <- u * t = u * x1 * x2
	fe_mul_u(r.t, r.t);

	// t2b <- t2b - y = (x1 + y1) (x2 + y2) - x1 * x2 - y1 * y2 = (x1 * y2 + y1 * x2)
	fe_sub(t2b, r.y, t2b);

	// z <- z1 * z2
	if (!z2_one) {
		fe_mul(p1.z, p2.z, r.z);
	}

	// t <- y + t = y1 * y2 + u * x1 * x2
	fe_add(r.y, r.t, r.t);

	// w1 <- z - w2 = z1 * z2 - d * u * t1 * t2
	if (z2_one) {
		fe_sub(p1.z, w2, w1);
	} else {
		fe_sub(r.z, w2, w1);
	}

	// z <- z + w2 = z1 * z2 + d * u * t1 * t2
	if (z2_one) {
		fe_add(p1.z, w2, r.z);
	} else {
		fe_add(r.z, w2, r.z);
	}

	// x3 <- t2b * w1 = (x1 * y2 + y1 * x2) * (z1 * z2 - d * u * t1 * t2)
	fe_mul(t2b, w1, r.x);

	// y3 <- z * t = (z1 * z2 + d * u * t1 * t2) * (y1 * y2 + u * x1 * x2)
	fe_mul(r.z, r.t, r.y);

	// t3 <- t2b * t = (x1 * y2 + y1 * x2) * (y1 * y2 + u * x1 * x2)
	if (out_precomp_t3) {
		fe_mul(r.t, t2b, r.t);
	}
	// Else: Not performed here: Can be reconstructed from t3 = r.t * t2b later

	// z3 <- w1 * z = (z1 * z2 - d * u * t1 * t2) * (z1 * z2 + d * u * t1 * t2)
	fe_mul(w1, r.z, r.z);
}

// Compute affine coordinates for (X, Y) from (X : Y : Z)
static void ec_affine(const ecpt &a, ecpt_affine &r) {
	// B = 1 / in.Z
	ufe b;
	fe_inv(a.z, b);

	// out.X = B * in.X
	fe_mul(a.x, b, r.x);

	// out.Y = B * in.Y
	fe_mul(a.y, b, r.y);

	// Final reduction
	fe_complete_reduce(r.x);
	fe_complete_reduce(r.y);
}

/*
 * Input validation:
 *
 * When the input point is not validated or other countermeasures are not
 * in place, it is possible to provide an input point on the twist of the
 * curve.  As shown in [7] this can lead to an active attack on the
 * cryptosystem.
 *
 * Bernstein's Curve25519 [8] prevents this attack by being "twist-secure",
 * for example, rather than validating the input.
 *
 * To avoid any invalid point fault attacks in my cryptosystem, I validate
 * that the input point (x, y) is on the curve, which is a cheap operation.
 *
 * I further check that the point is not x = 0, which would be another way
 * to introduce a fault, since x = 0 is the identity element.
 *
 * The input needs to fit within the field, so the exceptional value of
 * 2^127-1 must be checked for, since it is equivalent to 0.
 *
 * It is *not* verified that the y coordinate is non-zero.  When y = 0,
 * the point is in a cofactor subgroup rather than our q torsion group
 * of the base point.  Since these points can be introduced at any time
 * during the computation, checking for them would add unnecessary
 * complexity.  Instead, the ecmul functions double the output points
 * twice, preventing an attacker who injects invalid points of this
 * type from gaining any information about the rest of the scalar bits.
*/

// Verify that the affine point (x, y) exists on the given curve
// WARNING: Not constant time
static bool ec_valid_vartime(const ecpt_affine &p) {
	// If the platform is not little-endian,
	if (!IsLittleEndian()) {
		// We cannot work on this platform
		return false;
	}

	// If point is outside of field,
	if (!fe_infield_vartime(p.x) || !fe_infield_vartime(p.y)) {
		return false;
	}

	// If point is the additive identity x=0,
	if (fe_iszero_vartime(p.x)) {
		return false;
	}

	// Validate that:
	// 0 = u * x^2 - y^2 + 1 + d * u * x^2 * y^2

	// b <- x^2
	ufe b;
	fe_sqr(p.x, b);

	// b <- u * b = u * x^2
	fe_mul_u(b, b);

	// c <- y^2
	ufe c;
	fe_sqr(p.y, c);

	// Validate that:
	// 0 = b - c + 1 + d * b * c

	// r <- d * b
	ufe r;
	fe_mul_smallk(b, EC_D, r);

	// r <- r * c
	fe_mul(r, c, r);

	// r <- r + 1
	fe_add_smallk(r, 1, r);

	// r <- r - c
	fe_sub(r, c, r);

	// r <- r + b
	fe_add(r, b, r);

	// If the result is zero, it is on the curve
	return fe_iszero_vartime(r); // supports unreduced input
}

/*
 * Generate 64-bit mask:
 *
 * If x == y: return -1
 * Else: return 0
 */

// x =?= y
static CAT_INLINE s64 ec_gen_mask(const int x, const int y) {
	// Produce a 32-bit mask:
	// if x == y => -1
	// if x != y => 0
	s32 m32 = 0;
	m32 -= x ^ y;
	m32 >>= 31;
	m32 = ~m32;

	// Sign-extend it to 64 bits:
	return (s64)m32;
}

/*
 * Conditionally negate a point:
 *
 * If bit == 0: No change
 * Else if bit == 1: Negate the point
 *
 * Preconditions: bit is {0, 1}
 */

// r = [-]r
static CAT_INLINE void ec_cond_neg_inplace(const s32 bit, ecpt &r) {
	// Generate mask = -1 when bit == 1, else 0
	const u64 mask = (s64)(-bit);

	// Effectively negate when bit == 1
	ec_neg_mask_inplace(mask, r);
}

// r = [-]p
static CAT_INLINE void ec_cond_neg(const s32 bit, const ecpt &p, ecpt &r) {
	// Generate mask = -1 when bit == 1, else 0
	const u64 mask = (s64)(-bit);

	// Effectively negate when bit == 1
	ec_neg_mask(mask, p, r);
}

// r = [-]r
static CAT_INLINE void ec_cond_neg_affine(const s32 bit, ecpt_affine &r) {
	// Generate mask = -1 when bit == 1, else 0
	const u64 mask = (s64)(-bit);

	// Effectively negate when bit == 1
	ec_neg_mask_affine(mask, r);
}

/*
 * Conditionally add a point:
 *
 * If bit == 0: No change
 * Else if bit == 1: Returns sum of points
 *
 * Preconditions: bit is {0, 1}
 */
static CAT_INLINE void ec_cond_add(const s32 bit, const ecpt &a, const ecpt &b, ecpt &r, bool z2_one, bool precomp_t1, ufe &t2b) {
	// Generate mask = -1 when bit == 1, else 0
	const u64 mask = (s64)(-bit);

	// Set temporary point as identity or provided point
	ecpt T;
	ec_identity(T);
	ec_set_mask(b, mask, T);

	// Always add generated point
	ec_add(a, T, r, z2_one, precomp_t1, false, t2b);
}

/*
 * Elligator point decoding [18]
 *
 * This function "decodes" a random 256-bit integer into a deterministicly
 * chosen point on the curve, for which the private key is likely unknown.
 *
 * The input buffer is used as follows, taken as four 64-bit numbers N0..N3:
 *
 * a = (((N1 >> 1) << 64) | N0) + i * (((N3 >> 1) << 64) | N2).
 * sign_bit = Low bit of N3.
 *
 * Note that the low bit of N1 is entirely unused, and the low bit of N3 is
 * used to twiddle the sign of the curve point X coordinate.
 *
 * a*x^2 + y^2 = 1 + d*x^2*y^2
 * a = -u
 * d = 109*u
 * 
 * t^2 = s^3 + A*s^2 + B*s
 * A = -(a + d) = -(-1+109)u = -108u
 * B = ad = -109u^2
 *
 * z = -A / (1 + u*a^2)
 * e = chi(z^3 + A*z^2 + B*z)
 * s = e*z - (1-e)*A/2
 * t^2 = s^3 + A*s^2 + B*s
 *
 * The twisted Edwards curve Y coordinate is calculated first:
 *
 * Ey = (t^2 - (a-d)*s^2) / (t^2 + (a-d)*s^2)
 *
 * And then the curve equation is used to solve for the X coordinate:
 *
 * Ex = sqrt((Ey^2 - 1) / (d*Ey^2 - a))
 *
 * The sign of Ex is flipped based on the sign_bit described earlier.
 */
static void ec_elligator_decode(const char a0[32], ecpt_affine &r) {
	// Unpack random bytes into endian-neutral words
	ufe a;
	const u64 *words = reinterpret_cast<const u64 *>( a0 );
	a.a.i[0] = getLE(words[0]);
	a.a.i[1] = getLE(words[1]) >> 1;
	a.b.i[0] = getLE(words[2]);

	// Store final low bit of high word as a -1 or 0 mask
	u64 high_mask = getLE(words[3]);
	a.b.i[1] = high_mask >> 1;
	high_mask = -(s64)(high_mask & 1);

	// Note that one bit is ignored from the input.

	// There is a rare case where the random input will be 7FF..FFF,
	// which is the same as 0 in the base finite field Fp.
	// This rare case causes a slight acceptable bias towards (0,0).

	// Otherwise the input is uniformly distributed, and since
	// p = 1 mod 4, all input values in the field will map to curve
	// points according to [18].

	// It may be possible for the generated curve point to be invalid
	// as input for Snowshoe functions in rare cases, so it would be
	// a good idea to validate the output point before using it, and
	// generating new random input in this case.

	// z = -A / (1 + u * a^2)
	ufe z;
	fe_sqr(a, z);
	fe_mul_u(z, z);
	fe_add_smallk(z, 1, z);
	fe_inv(z, z);
	fe_mul_u(z, z);
	fe_mul_smallk(z, 108, z);

	// z2 = z^2
	ufe z2;
	fe_sqr(z, z2);

	// w1 = 108 * u * z^2
	ufe w1;
	fe_mul_u(z2, w1);
	fe_mul_smallk(w1, 108, w1);

	// z = 109 * u * z
	ufe w3;
	fe_mul_u(z, w3);
	fe_mul_u(w3, w3);
	fe_mul_smallk(w3, 109, w3);

	// e = z^3 - 108 * u * z^2 - 109 * u * u * z
	ufe e;
	fe_mul(z2, z, e);
	fe_sub(e, w1, e);
	fe_sub(e, w3, e);

	// e = chi(z^3 + A*z^2 + B*z)
	int chi = fe_chi(e);

	// s = e*z - (1 - e)*A/2 = z
	ufe s;

	// s = -z + 108u
	fe_set_smallk(108, s);
	fe_mul_u(s, s);
	fe_sub(s, z, s);

	// If chi == 1, s = z (constant-time)
	u64 mask = ~(s64)(chi >> 1);
	fe_set_mask(z, mask, s);

	// s2 = s^2
	ufe s2;
	fe_sqr(s, s2);

	// w2 = 108 * u * s^2
	ufe w2;
	fe_mul_u(s2, w2);
	fe_mul_smallk(w2, 108, w2);

	// s = 109 * u * u * s
	fe_mul_u(s, w3);
	fe_mul_u(w3, w3);
	fe_mul_smallk(w3, 109, w3);

	// t^2 = s^3 - 108 * u * s^2 - 109 * u * u * s
	ufe t2;
	fe_mul(s2, s, t2);
	fe_sub(t2, w2, t2);
	fe_sub(t2, w3, t2);

	// s2 = 110 * u * s^2
	fe_mul_u(s2, s2);
	fe_mul_smallk(s2, 110, s2);

	// d = 1 / (t^2 - 110 * u * s^2)
	ufe d;
	fe_sub(t2, s2, d);
	fe_inv(d, d);

	// r.y = (t^2 + 110 * u * s^2) / (t^2 - 110 * u * s^2)
	ufe y;
	fe_add(t2, s2, y);
	fe_mul(y, d, y);
	fe_set(y, r.y);

	// y2 = y^2
	ufe y2;
	fe_sqr(y, y2);

	// x = 1 / ((109 * y^2 + 1) * u)
	ufe x;
	fe_mul_smallk(y2, 109, x);
	fe_add_smallk(x, 1, x);
	fe_mul_u(x, x);
	fe_inv(x, x);

	// r.x = sqrt((y^2 - 1) / ((109 * y^2 + 1) * u))
	ufe f;
	fe_sub_smallk(y2, 1, f);
	fe_mul(f, x, x);
	fe_sqrt(x, x, false);

	// r.x = [-]x, based on one of the random input bits
	// It is clear from the visual shape of a twisted Edwards curve and the curve
	// equation that +X,-X and +Y,-Y are all valid points (hence group order cofactor 4).
	// Since the fe_sqrt() produces a somewhat unreliable sign, we are free to pick one
	// at random based on the input and both are valid.
	fe_neg_mask(high_mask, x, r.x);
}

