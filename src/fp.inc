#include "fp.hpp"
using namespace cat;

/*
 * 127-bit F(p) finite field arithmetic
 *
 * This is simply bigint math modulo Mersenne prime p = (2^127 - 1),
 * which admits perhaps the simplest, timing-invariant reduction.
 *
 * The inputs to the math functions are assumed to be partially-reduced,
 * meaning that I only assume that the high bit is clear, but they can
 * take on a value of 2^127-1 = p, which is equivalent to 0 in this
 * finite field.  The outputs of these functions will also be partially-
 * reduced.
 *
 * Before storing the results of these math functions, it is important to
 * call fp_complete_reduce() to take the partially reduced intermediate
 * results and produce completey reduced output less than p.
 */

// Load ufp from endian-neutral data bytes (16)
static void fp_load(const u8 *x, ufp &r) {
	r.i[0] = getLE64(*(u64*)x);
	r.i[1] = getLE64(*(u64*)(x + 8));
}

// Save ufp to endian-neutral data bytes (16)
static void fp_save(const ufp x, u8 *r) {
	*(u64*)r = getLE64(x.i[0]);
	*(u64*)(r + 8) = getLE64(x.i[1]);
}

// Check if a == b
// WARNING: Not constant-time
// WARNING: Inputs must be completely reduced
static CAT_INLINE bool fp_isequal_vartime(const ufp a, const ufp b) {
	return (a.i[0] == b.i[0]) && (a.i[1] == b.i[1]);
}

// Check if a == b in constant-time
// WARNING: Inputs must be completely reduced
static CAT_INLINE bool fp_isequal_ct(const ufp a, const ufp b) {
	u64 zero = a.i[0] ^ b.i[0];
	   zero |= a.i[1] ^ b.i[1];

	u32 z = (u32)zero | (u32)(zero >> 32);
	return z == 0;
}

// Check if x is zero; supports unreduced input
// WARNING: Not constant-time
static CAT_INLINE bool fp_iszero_vartime(const ufp x) {
	// If x == 2^127-1,
	if (x.i[0] == 0xffffffffffffffffULL &&
		x.i[1] == 0x7fffffffffffffffULL) {
		return true;
	}

	return (x.i[0] == 0) && (x.i[1] == 0);
}

// Check if x is zero
// WARNING: Input must be completely reduced
static CAT_INLINE bool fp_iszero_ct(const ufp x) {
	const u64 zero = x.i[0] | x.i[1];

	const u32 z = (u32)zero | (u32)(zero >> 32);
	return z == 0;
}

// Verify that 0 <= x < p
// WARNING: Not constant-time
static CAT_INLINE bool fp_infield_vartime(const ufp x) {
	// If high bit is set,
	if (x.i[1] >> 63) {
		// Not in field
		return false;
	}

	// If x == 2^127-1,
	if (x.i[0] == 0xffffffffffffffffULL &&
		x.i[1] == 0x7fffffffffffffffULL) {
		// Not in field
		return false;
	}

	return true;
}

// r = k
static CAT_INLINE void fp_set_smallk(const u32 k, ufp &r) {
	u128_set(r.w, k);
}

// r = 0
static CAT_INLINE void fp_zero(ufp &r) {
	u128_set(r.w, 0);
}

// r = a
static CAT_INLINE void fp_set(const ufp a, ufp &r) {
	r.w = a.w;
}

// r = (mask == -1) ? a : r
static CAT_INLINE void fp_set_mask(const ufp a, const u64 mask, ufp &r) {
	// Given a mask that is either -1 or 0,
	//	when the mask is zero, this function does nothing.
	//	when the mask is -1, it replaces `r` with `a`.
	r.i[0] ^= (a.i[0] ^ r.i[0]) & mask;
	r.i[1] ^= (a.i[1] ^ r.i[1]) & mask;
}

// r ^= a & mask
static CAT_INLINE void fp_xor_mask(const ufp &a, const u64 mask, ufp &r) {
	// Given a mask that is either -1 or 0,
	//	when the mask is zero, this function does nothing.
	//	when the mask is -1, it XORs `a` into `r`.
	r.i[0] ^= a.i[0] & mask;
	r.i[1] ^= a.i[1] & mask;
}

/*
 * Reduction after addition mod 2^127-1:
 *
 * When the high bit is set, we need to subtract the modulus to bring it back
 * under 2^128.
 *
 * Subtracting the modulus p = 2^127 - 1 is the same as turning off the high
 * bit and adding one.  Can this cause the high bit to be set again?  The
 * answer is no.  Proof:
 *
 * The inputs are assumed to be partially reduced, so the largest input
 * values is p, which is technically outside of Fp and equivalent to 0.
 * In this worst case, p + p - p = p, which is also partially reduced.
 *
 * Therefore, we just have to toggle the high bit off and add 1 if it was set.
 * This can be achieved with three Intel instructions!
 */

static CAT_INLINE void fp_add_reduce(ufp &x) {

#if defined(CAT_ASM_ATT) && defined(CAT_WORD_64) && defined(CAT_ISA_X86)

	CAT_ASM_BEGIN
		"btrq $63, %1\n\t"
		"adcq $0, %0\n\t"
		"adcq $0, %1"
		: "+r" (x.i[0]), "+r" (x.i[1])
		:
		: "cc"
	CAT_ASM_END

#else

	u32 msb = (u32)(u128_high(x.w) >> 63);
	u128_clear_msb(x.w);
	u128_add(x.w, msb);

#endif

}

/*
 * Reduction after subtraction mod 2^127-1:
 *
 * r = a - b
 *
 * Assume inputs are partially reduced, in the range [0, p] inclusive.
 *
 * Split this into cases to check validity of the code:
 *
 * POSITIVE CASE: a >= b
 *
 * The result does not need reduction in this case and the high bit is clear.
 *
 * NEGATIVE CASE: a < b
 *
 * The high bit will be set in this case since the high bit is clear in each
 * of the inputs, and a borrow will happen into the high bit position.
 *
 * In the negative case, we will want to add the modulus back in.  Adding the
 * modulus is the same as clearing the high bit and subtracting 1.  But can
 * another borrow occur as a result?  Answer: Thankfully no, and the proof:
 *
 * As 'a' gets larger, there will always be at least 1 bit set before the
 * high bit, so the subtraction of 1 for reduction will halt borrowing before
 * hitting the high bit.
 *
 * 0 - p is the worst case scenario here.  In this scenario, the result will
 * be a high bit and low bit set and all other bits will be cleared.
 * Clearing the resulting high bit and subtracting 1 will not cause another
 * borrow.  The result will be 0 as expected.
 *
 * Another example at the other extreme is 0 - 1.  In this example all the
 * resulting bits are set.  The high bit is set as expected, so reduction
 * logic will trigger to clear the high bit and subtract 1.  The result will
 * be a clear high bit, all bits set to 1 except the low bit.  This is equal
 * to p - 1, which is how to represent -1 in this finite field as expected.
 */

static CAT_INLINE void fp_sub_reduce(ufp &x) {

#if defined(CAT_ASM_ATT) && defined(CAT_WORD_64) && defined(CAT_ISA_X86)

	CAT_ASM_BEGIN
		"btrq $63, %1\n\t"
		"sbbq $0, %0\n\t"
		"sbbq $0, %1"
		: "+r" (x.i[0]), "+r" (x.i[1])
		:
		: "cc"
	CAT_ASM_END

#else

	u32 msb = (u32)(u128_high(x.w) >> 63);
	u128_clear_msb(x.w);
	u128_sub(x.w, msb);

#endif

}

/*
 * Negation:
 *
 * This is the same as subtraction with r = 0 - a.
 * And reduction is performed the same way.
 *
 * Note that fp_neg(p) == 0 and fp_neg(0) == 0
 */

// r = -a
static CAT_INLINE void fp_neg(const ufp a, ufp &r) {
	// Uses 1a 1r
	r.w = u128_neg(a.w);
	fp_sub_reduce(r);
}

// r = (mask==-1 ? r : -r)
static CAT_INLINE void fp_neg_mask_inplace(const u64 mask, ufp &r) {
	// s = -p
	ufp s;
	s.w = u128_neg(r.w);

#if defined(CAT_ASM_ATT) && defined(CAT_WORD_64) && defined(CAT_ISA_X86)

	const u32 m = static_cast<const u32>( mask );

	CAT_ASM_BEGIN
		"testl %2, %2\n\t"
		"cmovnzq %3, %0\n\t"
		"cmovnzq %4, %1"
		: "+r" (r.i[0]), "+r" (r.i[1])
		: "r" (m), "r" (s.i[0]), "r" (s.i[1])
		: "cc"
	CAT_ASM_END

#else

	// If mask is set, take 's', else take 'r'
	r.i[0] ^= (s.i[0] ^ r.i[0]) & mask;
	r.i[1] ^= (s.i[1] ^ r.i[1]) & mask;

#endif

	// Note that fp_sub_reduce does nothing for the case of mask==0, since the
	// high bit will not be set in this case.
	fp_sub_reduce(r);
}

// r = (mask==-1 ? p : -p)
static CAT_INLINE void fp_neg_mask(const u64 mask, const ufp p, ufp &r) {
	// Assuming mask is either -1 or 0,
	//	when mask is zero, it sets r = p
	//	when mask is -1, it sets r = ~p + 1 = -p
	r.i[0] = p.i[0] ^ mask;
	r.i[1] = p.i[1] ^ mask;
	u128_add(r.w, mask & 1);

	// Note that fp_sub_reduce does nothing for the case of mask==0, since the
	// high bit will not be set in this case.
	fp_sub_reduce(r);
}

/*
 * Reduction from p to 0 is slightly sped up by negating twice.
 *
 * Note that fp_neg(p) => 0 and fp_neg(0) => 0.
 *
 * This takes ops: 8 adds, 2 btrq 
 *
 * Another way to implement it is by ANDing together all the bits
 * and then adding the resulting 1/0 bit to nudge p to 0.  This
 * takes 8 ands, 3/4 shifts, 4 adds, 1 btrq -- quite a bit slower.
 */

// Reduce r < p in the case where r = p
static CAT_INLINE void fp_complete_reduce(ufp &r) {
	// Uses 2a 2r
	fp_neg(r, r);
	fp_neg(r, r);
}

// r = a + b
static CAT_INLINE void fp_add(const ufp &a, const ufp &b, ufp &r) {
	// Uses 1a 1r
	r.w = u128_sum(a.w, b.w);
	fp_add_reduce(r);
}

// r = a - b
static CAT_INLINE void fp_sub(const ufp &a, const ufp &b, ufp &r) {
	// Uses 1a 1r
	r.w = u128_diff(a.w, b.w);
	fp_sub_reduce(r);
}

// r = a + (u32)k
static CAT_INLINE void fp_add_smallk(const ufp &a, const u32 k, ufp &r) {
	// Uses 1a 1r
#ifdef CAT_HAS_U128
	r.w = a.w + k;
#else
	r.w = u128_sum(a.w, (u64)k);
#endif
	fp_add_reduce(r);
}

// r = a - (u32)k
static CAT_INLINE void fp_sub_smallk(const ufp &a, const u32 k, ufp &r) {
	// Uses 1a 1r
#ifdef CAT_HAS_U128
	r.w = a.w - k;
#else
	r.w = u128_diff(a.w, (u64)k);
#endif
	fp_sub_reduce(r);
}

/*
 * Multiplication of two partially-reduced inputs:
 *
 * The annotations in the comments prove at each step that any partially-reduced
 * input will not overflow any of the intermediate results.
 *
 * Using the schoolbook method:
 *
 *             B1 B0
 *           x A1 A0
 *          --------
 *             00 00 <- low product     = A0*B0
 *          01 01    <- middle products = A1*B0 + A0*B1
 *     + 11 11       <- high product    = A1*B1
 *     -------------
 *
 * The middle products' sum fits in 128 bits so I can skip a reduction step.
 *
 * The high half of the low and middle products can be rolled into the next
 * product up without overflow, so the low parts of each are left behind.
 *
 * Then the high product is shifted left, the high bit of the low half of
 * the 256-bit product is stolen into the high part, and then the high half
 * of the 256-bit product is subtracted from the low half to reduce it to
 * 127 bits.
 *
 * This reduction approach takes full advantage of the partially reduced
 * inputs to avoid extra operations and appears to be minimal.
 */

// r = a * b
static CAT_INLINE void fp_mul(const ufp a, const ufp b, ufp &r) {
	// Uses 4m 5a 1r
	// a.i[0] = A0, a.i[1] = A1, b.i[0] = B0, b.i[1] = B1

	// middle = A0*B1 + B1*B0
	// middle <= 2*(2^64-1)(2^63-1)
	//         = (2^64-1)(2^64-2)
	//         = 2^128 - 3 * 2^64 + 2
	//        <= 2^128 - 1
	// Sum fits within 128 bits
	u128 middle = u128_prod(u128_low(a.w), u128_high(b.w));
	u128_add(middle, u128_prod(u128_low(b.w), u128_high(a.w)));

	// low = A0*B0 <= (2^128-1)
	u128 low = u128_prod(u128_low(a.w), u128_low(b.w));

	// Incorporate high half of low product into the middle
	// product, so that it can be neglected:

	// middle += high half of low
	// middle <= (2^64-1) + (2^128 - 3 * 2^64 + 2)
	//         = 2^128 - 2 * 2^64 + 1
	//        <= 2^128 - 1
	// Sum fits within 128 bits
	u128_add(middle, u128_high(low));

	// high = A1*B1
	// high <= (2^63-1)(2^63-1)
	//       = 2^126 - 2 * 2^63 + 1
	//       = 2^126 - 2^64 + 1
	// Product is completely reduced
	ufp high;

	// Incorporate the high half of middle product into the
	// high product, so that it can be neglected:

	// high += high half of middle
	// high <= 2^126 - 2^64 + 1 + (2^64 - 2)
	//       = 2^126 - 1
	// Sum is completely reduced

	high.w = u128_prod_63_sum(u128_high(a.w), u128_high(b.w), u128_high(middle));

	/*
	 * The situation is now:
	 *
	 *             B1 B0
	 *           x A1 A0
	 *          --------
	 *                00 <- Only low part remains (incorporated above)
	 *             01    <- Only low part remains (incorporated above)
	 *     + 11 11       <- high terms   = A1*B1
	 *     -------------
	 */

	// To reduce the 256-bit product, we need to multiply the binary number
	// represented by all of the overflow bits from 127 and up by 2^127-1,
	// and then subtract it from the low 128 bits of the product.  The high
	// product can be doubled, which still fits within 128 bits.  The high
	// bit of the low 128 bits can be added to it to create the number
	// represented by the overflow bits.

	// Incorporate high bit of middle into shifted high bits:

	// high = 2 * high + (high bit of low part of the middle product)
	// high <= 2^127 - 2 + (1)
	//       = 2^127 - 1
	// Sum is partially reduced
	u128_lshift(high.w, 1);
	u128_or(high.w, u128_low(middle) >> 63);

	// Store low 127 bits in result, as high bit was incorporated above.
	// Result is partially reduced
	u128_set(r.w, u128_low(low), u128_low(middle) & 0x7fffffffffffffffULL);

	// Subtract p * high from r:
	// This is the same as adding high to r, where both
	// high and r happen to be partially reduced inputs.

	// Add two partially reduced inputs with proven Fp addition:
	fp_add(high, r, r);
}

// r = a * b, b = small 32-bit constant
static CAT_INLINE void fp_mul_smallk(const ufp a, const u32 b, ufp &r) {
	// Uses 2m 3a 1r

	// Eliminate multiplications by high part of b, which are 0 in this
	// special case:

	u128 low = u128_prod(u128_low(a.w), b);
	// low <= (2^64-1)(2^32-1)
	//      = 2^96 - 2^64 - 2^32 + 1

	// middle = (high half of a) * b
	// middle <= (2^63-1)(2^32-1)
	//         = 2^95 - 2^63 - 2^32 + 1

	// middle += (high half of low)
	// middle <= 2^95 - 2^63 - 2^32 + 1 - (2^32) - (1)
	//         = 2^95 - 2^63 - 2^33

	u128 middle = u128_prod_63_sum(u128_high(a.w), b, u128_high(low));

	// Only calculate a 64-bit sum for high in this special case:

	ufp high;
	u128_set(high.w, (u128_high(middle) << 1) | (u128_low(middle) >> 63));
	// high <= 2^32

	u128_set(r.w, u128_low(low), u128_low(middle) & 0x7fffffffffffffffULL);

	fp_add(high, r, r);
}

// r = a^2
static CAT_INLINE void fp_sqr(const ufp a, ufp &r) {
	// Uses 3m 5a 1r

	// In this special case the cross terms are equal, so
	// trade one multiplication for a shift:

	const u64 x = u128_low(a.w);
	const u64 y = u128_high(a.w);

	u128 middle = u128_prod(x, y);
	u128_lshift(middle, 1);

	// The rest is the same as fp_mul:

	u128 low = u128_prod(x, x);

	u128_add(middle, u128_high(low));

	ufp high;
	high.w = u128_prod_63_sum(y, y, u128_high(middle));

	u128_lshift(high.w, 1);
	u128_or(high.w, u128_low(middle) >> 63);

	u128_set(r.w, u128_low(low), u128_low(middle) & 0x7fffffffffffffffULL);

	fp_add(high, r, r);
}

// r = a / 2
static CAT_INLINE void fp_div2(const ufp a, ufp &r) {
	// Uses 5a 1r

	// This is derived from fp_mul, with 0x4000...00,
	// which is (1/2) in Fp.

	const u64 x = u128_low(a.w);
	const u64 y = u128_high(a.w);

	u128 middle;
	u128_set(middle, x << 62, x >> 2);

	ufp high;
	u128_set(high.w, y << 62, y >> 2);

	u128_add(high.w, u128_high(middle));

	u128_lshift(high.w, 1);
	u128_or(high.w, u128_low(middle) >> 63);

	u128_set(r.w, 0, u128_low(middle) & 0x7fffffffffffffffULL);

	fp_add(high, r, r);
}

// r = 1/x
static void fp_inv(const ufp x, ufp &r) {
	// Uses 126S 12M

	/*
	 * Euler's totient function:
	 * 1/x = x ^ (2^127 - 1 - 2)
	 *
	 * Evaluating a potential alternative: Use a short addition chain?
	 *
	 * l(127): 1 2 3 6 12 24 48 51 63 126 127 = 10 ops for 7 bits
	 * l(2^n-1) ~= n + l(n) - 1, n = 127 + 10 - 1 = 136
	 *
	 * Straight-forward squaring takes 138 ops.  Seems just fine.
	 */
	ufp n1, n2, n3, n4, n5, n6;

	fp_sqr(x, n2);
	fp_mul(x, n2, n2);
	fp_sqr(n2, n3);
	fp_sqr(n3, n3);
	fp_mul(n3, n2, n3);
	fp_sqr(n3, n4); fp_sqr(n4, n4);
	fp_sqr(n4, n4); fp_sqr(n4, n4);
	fp_mul(n3, n4, n4);
	fp_sqr(n4, n5); fp_sqr(n5, n5); fp_sqr(n5, n5); fp_sqr(n5, n5);
	fp_sqr(n5, n5); fp_sqr(n5, n5); fp_sqr(n5, n5); fp_sqr(n5, n5);
	fp_mul(n5, n4, n5);
	fp_sqr(n5, n6); fp_sqr(n6, n6); fp_sqr(n6, n6); fp_sqr(n6, n6);
	fp_sqr(n6, n6); fp_sqr(n6, n6); fp_sqr(n6, n6); fp_sqr(n6, n6);
	fp_sqr(n6, n6); fp_sqr(n6, n6); fp_sqr(n6, n6); fp_sqr(n6, n6);
	fp_sqr(n6, n6); fp_sqr(n6, n6); fp_sqr(n6, n6); fp_sqr(n6, n6);
	fp_mul(n5, n6, n6);
	fp_sqr(n6, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_mul(n1, n6, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_mul(n1, n6, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_mul(n1, n5, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_mul(n1, n4, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_mul(n1, n3, n1);
	fp_sqr(n1, n1);
	fp_mul(n1, x, n1);
	fp_sqr(n1, n1);
	fp_sqr(n1, n1);
	fp_mul(n1, x, r);
}

// r = sqrt(x)
static void fp_sqrt(const ufp x, ufp &r) {
	// Uses 125S

	// sqrt(x) = x ^ ((p+1)/4) = x ^ (2^125)
	ufp t;
	fp_sqr(x, t);
	for (int ii = 1; ii < 124; ++ii) {
		fp_sqr(t, t);
	}
	fp_sqr(t, r);
}

// r = chi(x)
// -1 if 'x' does not have a square root.
//  0 is 'x' is zero.
// +1 if 'x' has a square root.
static int fp_chi(const ufp x) {
	// Uses 126S 11M

	/*
	 * chi(x) = x ^ ((p-1)/2) = x ^ (2^126 - 1)
	 *
	 * This could also be done with the GCD but then it would not
	 * be constant-time.
	 */
	ufp n1, n2, n3, n4, n5, n6;

	fp_sqr(x, n2);
	fp_mul(x, n2, n2);
	fp_sqr(n2, n3);
	fp_sqr(n3, n3);
	fp_mul(n3, n2, n3);
	fp_sqr(n3, n4); fp_sqr(n4, n4);
	fp_sqr(n4, n4); fp_sqr(n4, n4);
	fp_mul(n3, n4, n4);
	fp_sqr(n4, n5); fp_sqr(n5, n5); fp_sqr(n5, n5); fp_sqr(n5, n5);
	fp_sqr(n5, n5); fp_sqr(n5, n5); fp_sqr(n5, n5); fp_sqr(n5, n5);
	fp_mul(n5, n4, n5);
	fp_sqr(n5, n6); fp_sqr(n6, n6); fp_sqr(n6, n6); fp_sqr(n6, n6);
	fp_sqr(n6, n6); fp_sqr(n6, n6); fp_sqr(n6, n6); fp_sqr(n6, n6);
	fp_sqr(n6, n6); fp_sqr(n6, n6); fp_sqr(n6, n6); fp_sqr(n6, n6);
	fp_sqr(n6, n6); fp_sqr(n6, n6); fp_sqr(n6, n6); fp_sqr(n6, n6);
	fp_mul(n5, n6, n6);
	fp_sqr(n6, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_mul(n1, n6, n1); // n1 = 2^64 - 1
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_mul(n1, n6, n1); // n1 = 2^96 - 1
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_mul(n1, n5, n1); // n1 = 2^112 - 1
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_mul(n1, n4, n1); // n1 = 2^120 - 1
	fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1); fp_sqr(n1, n1);
	fp_mul(n1, n3, n1); // n1 = 2^124 - 1
	fp_sqr(n1, n1);
	fp_sqr(n1, n1);
	fp_mul(n1, n2, n1); // n1 = 2^126 - 1

	// Note that now n1 =
	// "-1" is 7fff.fffe
	// "0" is 0
	// "1" is 1

	const u64 n = u128_low(n1.w);

	// Convert to standard integer -1
	return (s32)(n >> 1) | (s32)n;
}

