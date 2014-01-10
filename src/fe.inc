// Based on Fp finite field
#include "fp.inc"
#include "fe.hpp"

/*
 * 254-bit GF(p^2) optimal extension field (OEF) arithmetic
 *
 * This is simply bigint complex math (a + ib), with a,b in Fp.
 *
 * Before storing the results of these math functions, it is important to
 * call fe_complete_reduce() to produce fully-reduced output.
 */

// Load ufe from endian-neutral data bytes (32)
static void fe_load(const u8 *x, ufe &r) {
	fp_load(x, r.a);
	fp_load(x + 16, r.b);
}

// Save ufe to endian-neutral data bytes (32)
static void fe_save(const ufe &x, u8 *r) {
	fp_save(x.a, r);
	fp_save(x.b, r + 16);
}

// Check if r is zero; supports unreduced input
// WARNING: Not constant-time
static CAT_INLINE bool fe_iszero_vartime(const ufe &r) {
	return fp_iszero_vartime(r.a) && fp_iszero_vartime(r.b);
}

// Check if r is zero in constant-time
// WARNING: Input must be completely reduced
static CAT_INLINE bool fe_iszero_ct(const ufe &r) {
	u64 zero = r.a.i[0] | r.a.i[1] | r.b.i[0] | r.b.i[1];
	u32 z = (u32)(zero | (zero >> 32));
	return z == 0;
}

// Check if a == b
// WARNING: Not constant-time
// WARNING: Inputs must be completely reduced
static CAT_INLINE bool fe_isequal_vartime(const ufe &a, const ufe &b) {
	return fp_isequal_vartime(a.a, b.a) && fp_isequal_vartime(a.b, b.b);
}

// Check if a == b
// WARNING: Inputs must be completely reduced
static CAT_INLINE bool fe_isequal_ct(const ufe &a, const ufe &b) {
	u64 zero = (a.a.i[0] ^ b.a.i[0]) | (a.a.i[1] ^ b.a.i[1]) |
			   (a.b.i[0] ^ b.b.i[0]) | (a.b.i[1] ^ b.b.i[1]);
	u32 z = (u32)zero | (u32)(zero >> 32);
	return z == 0;
}

// Validate that r is within field
// WARNING: Not constant-time
static CAT_INLINE bool fe_infield_vartime(const ufe &r) {
	return fp_infield_vartime(r.a) && fp_infield_vartime(r.b);
}

// Reduce r
static CAT_INLINE void fe_complete_reduce(ufe &r) {
	// If a or b = 2^127-1, set that one to 0.
	// NOTE: The math functions already ensure that a,b < 2^127

	fp_complete_reduce(r.a);
	fp_complete_reduce(r.b);
}

// r = (k + 0i)
static CAT_INLINE void fe_set_smallk(const u32 k, ufe &r) {
	fp_set_smallk(k, r.a);
	fp_zero(r.b);
}

// r = 0
static CAT_INLINE void fe_zero(ufe &r) {
	fp_zero(r.a);
	fp_zero(r.b);
}

// r = a
static CAT_INLINE void fe_set(const ufe &a, ufe &r) {
	fp_set(a.a, r.a);
	fp_set(a.b, r.b);
}

// r = (mask == -1) ? a : r
static CAT_INLINE void fe_set_mask(const ufe &a, const u64 mask, ufe &r) {
	fp_set_mask(a.a, mask, r.a);
	fp_set_mask(a.b, mask, r.b);
}

// r ^= a & mask
static CAT_INLINE void fe_xor_mask(const ufe &a, const u64 mask, ufe &r) {
	// I tried inline assembly code here for CMOV, which seems like the best
	// place to do it.  That unexpectedly reduced performance noticeably.
	fp_xor_mask(a.a, mask, r.a);
	fp_xor_mask(a.b, mask, r.b);
}

// r = a'
static CAT_INLINE void fe_conj(const ufe &a, ufe &r) {
	// Uses 1A

	fp_set(a.a, r.a);
	fp_neg(a.b, r.b);
}

// r = -a
static CAT_INLINE void fe_neg(const ufe &a, ufe &r) {
	// Uses 2A

	fp_neg(a.a, r.a);
	fp_neg(a.b, r.b);
}

// r = (mask==-1 ? r : -r)
static CAT_INLINE void fe_neg_mask_inplace(const u64 mask, ufe &r) {
	fp_neg_mask_inplace(mask, r.a);
	fp_neg_mask_inplace(mask, r.b);
}

// r = (mask==-1 ? p : -p)
static CAT_INLINE void fe_neg_mask(const u64 mask, const ufe &p, ufe &r) {
	fp_neg_mask(mask, p.a, r.a);
	fp_neg_mask(mask, p.b, r.b);
}

// r = r + (k + 0i)
static CAT_INLINE void fe_add_smallk(const ufe &a, const u32 k, ufe &r) {
	// Uses 1A

	fp_add_smallk(a.a, k, r.a);
	fp_set(a.b, r.b);
}

// r = a + b
static CAT_INLINE void fe_add(const ufe &a, const ufe &b, ufe &r) {
	// Uses 2A

	// Seems about comparable to 2^^256-c in performance
	fp_add(a.a, b.a, r.a);
	fp_add(a.b, b.b, r.b);
}

// r = a - b
static CAT_INLINE void fe_sub(const ufe &a, const ufe &b, ufe &r) {
	// Uses 2A

	// Seems about comparable to 2^^256-c in performance
	fp_sub(a.a, b.a, r.a);
	fp_sub(a.b, b.b, r.b);
}

// r = r - (k + 0i)
static CAT_INLINE void fe_sub_smallk(const ufe &a, const u32 k, ufe &r) {
	// Uses 1A

	fp_sub_smallk(a.a, k, r.a);
	fp_set(a.b, r.b);
}

// r = a * u, u = 2 + i
static CAT_INLINE void fe_mul_u(const ufe &a, ufe &r) {
	// (a0 + ia1) * (2 + i)
	// = (a0*2 - a1) + i(a1*2 + a0)

	// t0 <- a - b
	ufp t0;
	fp_sub(a.a, a.b, t0);

	// t1 <- b + a
	ufp t1;
	fp_add(a.b, a.a, t1);

	// a' <- a + a - b = t0 + a
	fp_add(t0, a.a, r.a);

	// b' <- b + b + a = t1 + b
	fp_add(t1, a.b, r.b);
}

// r = a * b
static void fe_mul(const ufe &a, const ufe &b, ufe &r) {
	// Uses 3M 5A

	// (a0 + ia1) * (b0 + ib1)
	// = (a0b0 - a1b1) + i(a1b0 + a0b1)
	// = (a0b0 - a1b1) + i( (a1 + a0) * (b1 + b0) - a1b1 - a0b0 )

	ufp t0, t1, t2, t3;

	// Ordered to take advantage of pipeline ILP:

	fp_add(a.a, a.b, t0);
	fp_add(b.a, b.b, t1);
	fp_mul(a.a, b.a, t2);
	fp_mul(t0, t1, t1);
	fp_mul(a.b, b.b, t3);
	fp_sub(t1, t2, t1);
	fp_sub(t2, t3, r.a);
	fp_sub(t1, t3, r.b);
}

// r = a * b(small constant)
static CAT_INLINE void fe_mul_smallk(const ufe &a, const u32 b, ufe &r) {
	// Uses 2m

	fp_mul_smallk(a.a, b, r.a);
	fp_mul_smallk(a.b, b, r.b);
}

// r = a ^ 2
static CAT_INLINE void fe_sqr(const ufe &a, ufe &r) {
	// Uses 2M 3A

	// (a + ib) * (a + ib)
	// = (aa - bb) + i(ab + ab)
	// = (a + b) * (a - b) + i(ab + ab)

	ufp t0, t1, t2;

	fp_add(a.a, a.b, t0);
	fp_sub(a.a, a.b, t1);
	fp_add(a.a, a.a, t2);
	fp_mul(t0, t1, r.a);
	fp_mul(t2, a.b, r.b);
}

// r = 1 / x
static void fe_inv(const ufe &x, ufe &r) {
	// Uses 2S 2M 2A 1FpInv

	// 1/x = x'/|x|
	// NOTE: The inversion only needs to be done over a 2^^127 field instead of 2^^256

	ufp t0, t1, t2;

	fp_sqr(x.a, t0);
	fp_sqr(x.b, t1);
	fp_add(t0, t1, t2);

	fp_inv(t2, t0);

	fp_neg(x.b, t1);

	fp_mul(x.a, t0, r.a);
	fp_mul(t1, t0, r.b);
}

// r = chi(x)
static int fe_chi(const ufe &x) {
	// Uses 2S 1A 1FpChi

	// chi(x) = x ^ ((p^2-1)/2)
	//        = x ^ ( (p-1)/2 * (1+p) )
	//        = |x| ^ ((p-1)/2)
	// This uses the identity x^(1+p) = |x|.
	// NOTE: The heavy lifting is done over a 2^127 field instead of 2^256

	// t0 = |x|
	ufp t0, t1;
	fp_sqr(x.a, t0);
	fp_sqr(x.b, t1);
	fp_add(t0, t1, t0);

	// chi(t0)
	return fp_chi(t0);
}

// r = sqrt(x)
// Note that the sign on the result is not necessarily sgn(x)
static bool fe_sqrt(const ufe &x, ufe &r, bool check_input_vartime) {
	// Requires 2FpSqrt 1FpInv 1FpChi, in constant-time

	// Uses a well-known algorithm, which is well stated in
	// "Square root computation over even extension fields" (Adj Henriquez 2012)
	// http://eprint.iacr.org/2012/685.pdf
	// They call it the "Complex method"; it's the fastest one in the paper.

	// Yes, Fp^2 is fucking terrible for square roots.
	// Fortunately it is only needed for Elligator.

	ufp alpha, delta, delta2, t;

	// Note that most of these operations can be skipped if
	// x.b == 0.  However that would make this variable-time.

	// alpha = x.a^2 + x.b^2
	fp_sqr(x.a, alpha);
	fp_sqr(x.b, t);
	fp_add(alpha, t, alpha);

	// If validating input,
	if (check_input_vartime) {
		// chi = chi(alpha) [expensive!]
		if (fp_chi(alpha) == -1) {
			return false;
		}
	}
	// Otherwise this is constant-time.

	// alpha = sqrt(alpha)
	fp_sqrt(alpha, alpha);

	// delta = (x.a + alpha) / 2
	fp_add(x.a, alpha, delta);
	fp_div2(delta, delta);

	// chi = chi(delta)
	int chi = fp_chi(delta);

	// delta2 = (x.a - alpha) / 2
	fp_sub(x.a, alpha, delta2);
	fp_div2(delta2, delta2);

	// If chi == -1, use delta2 instead
	u64 mask = (s64)(chi >> 1);
	fp_set_mask(delta2, mask, delta);

	// r.a = sqrt(delta)
	fp_sqrt(delta, r.a);

	// r.b = x.b / (2*r.a)
	fp_add(r.a, r.a, t);
	fp_inv(t, t);
	fp_mul(x.b, t, r.b);

	return true;
}

