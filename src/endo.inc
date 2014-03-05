// Optimal extension field Fp^2
#include "fe.inc"
#include "ecpt.hpp"

/*
 * GLS Endomorphism [12]
 *
 * In affine twisted Edwards coordinates:
 * endomorphism(x, y) = (wx * conj(x), conj(y)) = lambda * (x, y)
 */

// Derivation in README.md
static const ufe ENDO_WX = {
	{{	// Real part (a):
		0x695AB4D883DE0B89ULL,
		0x59F30C694ED33218ULL
	}},
	{{	// Imaginary part (b):
		0xD2B569B107BC1713ULL,
		0x33E618D29DA66430ULL
	}}
};

/*
 * Scalar decomposition modulo q:
 *
 * This algorithm uses no conditional branches, memory allocation,
 * or potentially timing-variant opcodes (like division).
 *
 * Constants:
 * A := 2^126 - 1; // Replace multiplication by A with shift and subtract
 * B := 0x62D2CF00A287A526; // Relatively small 64-bit constant
 * qround := q div 2 + 1;
 * qround = 0x7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFD3130A0A606E43E9E74DB471D84F00D3;
 *
 * a1 := A * k;
 * a2 := B * k;
 *
 * These divisions are implemented using multiplication using the technique
 * from [2], which is somewhat slow but allows for constant-time operation:
 * z1 := (a1 + qround) div q;
 * z2 := (a2 + qround) div q;
 *
 * The final step can have negative intermediate results after subtractions:
 * k1 := k - (z1 * A + z2 * B);
 * k2 := z1 * B - z2 * A;
 *
 * k = (k1 + k2 * l) mod q, ||k1||, ||k2>|| <= 126 bits
 */

// k = k1 + k2 * l
static void gls_decompose(const u64 k[4], s32 &k1sign, ufp &k1, s32 &k2sign, ufp &k2) {
	static const u64 B = 0x62D2CF00A287A526ULL;

	// k < q < 2^252

	static const u64 QROUND[4] = {
		0xE74DB471D84F00D3ULL,
		0xD3130A0A606E43E9ULL,
		0xFFFFFFFFFFFFFFFFULL,
		0x07FFFFFFFFFFFFFFULL
	};

	u64 a1[6], a2[5], t[5];
	u128 sum, prod;

	// a1 <- (k << 126) - k = A * k < 2^(252+126 = 378)
	sum = u128_diff((u64)0, k[0]);
	a1[0] = u128_low(sum);
	u128_borrow_add_sub(sum, k[0] << 62, k[1]);
	a1[1] = u128_low(sum);
	u128_borrow_add_sub(sum, (k[1] << 62) | (k[0] >> 2), k[2]);
	a1[2] = u128_low(sum);
	u128_borrow_add_sub(sum, (k[2] << 62) | (k[1] >> 2), k[3]);
	a1[3] = u128_low(sum);
	u128_borrow_add(sum, (k[3] << 62) | (k[2] >> 2));
	a1[4] = u128_low(sum);
	a1[5] = u128_high(sum) + (k[3] >> 2);

	// a2 <- B * k < 2^(252+63 = 315)
	prod = u128_prod(B, k[0]);
	a2[0] = u128_low(prod);
	prod = u128_prod_sum(B, k[1], u128_high(prod));
	a2[1] = u128_low(prod);
	prod = u128_prod_sum(B, k[2], u128_high(prod));
	a2[2] = u128_low(prod);
	prod = u128_prod_sum(B, k[3], u128_high(prod));
	a2[3] = u128_low(prod);
	a2[4] = u128_high(prod);

	// a1 += qround < 2^379
	sum = u128_sum(QROUND[0], a1[0]);
	a1[0] = u128_low(sum);
	u128_carry_add(sum, a1[1], QROUND[1]);
	a1[1] = u128_low(sum);
	u128_carry_add(sum, a1[2], QROUND[2]);
	a1[2] = u128_low(sum);
	u128_carry_add(sum, a1[3], QROUND[3]);
	a1[3] = u128_low(sum);
	u128_carry_add(sum, a1[4]);
	a1[4] = u128_low(sum);
	a1[5] = u128_high(sum) + a1[5];

	// a2 += qround < 2^316
	sum = u128_sum(QROUND[0], a2[0]);
	a2[0] = u128_low(sum);
	u128_carry_add(sum, a2[1], QROUND[1]);
	a2[1] = u128_low(sum);
	u128_carry_add(sum, a2[2], QROUND[2]);
	a2[2] = u128_low(sum);
	u128_carry_add(sum, a2[3], QROUND[3]);
	a2[3] = u128_low(sum);
	a2[4] = u128_high(sum) + a2[4];

	/*
	 * Using the Unsigned Division algorithm from section 4 of [2]:
	 *
	 * d = q < 2^252, l = 252
	 * n = a1 < 2^379, N = 379
	 *
	 * m' = 2^N * (2^l - d) / d + 1
	 * = 2^(N+l)/d - 2^N + 1
	 *
	 * t := (mp * n) div (2^N);
	 *
	 * s := t + ((n - t) div 2);
	 *
	 * quot := s div (2^(l-1));
	 *
	 * See magma_unsigned_division.txt for more details.
	 */

	// m' = 0x2CECF5F59F91BC1618B24B8E27B0FF2E7C49FC5B010B8E4474AF5BD870D4C42C
	static const u64 M1[4] = {
		0x74AF5BD870D4C42CULL,
		0x7C49FC5B010B8E44ULL,
		0x18B24B8E27B0FF2EULL,
		0x2CECF5F59F91BC16ULL
	};

	// t <- m' * a1 >> (379 = 64 * 5 + 59)

	// Comba multiplication: Right to left schoolbook column approach
#define MUL_DIGIT_START(ii, jj) \
	prod = u128_prod_sum(M1[ii], a1[jj], u128_low(sum)); \
	sum = u128_sum(u128_high(sum), u128_high(prod));

#define MUL_DIGIT(ii, jj) \
	prod = u128_prod_sum(M1[ii], a1[jj], u128_low(prod)); \
	u128_add(sum, u128_high(prod));

	prod = u128_prod(M1[0], a1[0]);

	prod = u128_prod_sum(M1[1], a1[0], u128_high(prod));
	u128_set(sum, u128_high(prod));
	MUL_DIGIT(0, 1);

	MUL_DIGIT_START(2, 0);
	MUL_DIGIT(1, 1);
	MUL_DIGIT(0, 2);

	MUL_DIGIT_START(3, 0);
	MUL_DIGIT(2, 1);
	MUL_DIGIT(1, 2);
	MUL_DIGIT(0, 3);

	MUL_DIGIT_START(3, 1);
	MUL_DIGIT(2, 2);
	MUL_DIGIT(1, 3);
	MUL_DIGIT(0, 4);

	MUL_DIGIT_START(3, 2);
	MUL_DIGIT(2, 3);
	MUL_DIGIT(1, 4);
	MUL_DIGIT(0, 5);

	t[0] = u128_low(prod);

	MUL_DIGIT_START(3, 3);
	MUL_DIGIT(2, 4);
	MUL_DIGIT(1, 5);

	t[1] = u128_low(prod);

	MUL_DIGIT_START(3, 4);
	MUL_DIGIT(2, 5);

	t[2] = u128_low(prod);

	prod = u128_prod_sum(M1[3], a1[5], u128_low(sum));

	t[3] = u128_low(prod);
	t[4] = u128_high(prod) + u128_high(sum);

#undef MUL_DIGIT_START
#undef MUL_DIGIT

	// t >>= 59
	t[0] = (t[0] >> 59) | (t[1] << 5);
	t[1] = (t[1] >> 59) | (t[2] << 5);
	t[2] = (t[2] >> 59) | (t[3] << 5);
	t[3] = (t[3] >> 59) | (t[4] << 5);
	t[4] >>= 59;

	// a1 -= t
	sum = u128_diff(a1[0], t[0]);
	a1[0] = u128_low(sum);
	u128_borrow_add_sub(sum, a1[1], t[1]);
	a1[1] = u128_low(sum);
	u128_borrow_add_sub(sum, a1[2], t[2]);
	a1[2] = u128_low(sum);
	u128_borrow_add_sub(sum, a1[3], t[3]);
	a1[3] = u128_low(sum);
	u128_borrow_add_sub(sum, a1[4], t[4]);
	a1[4] = u128_low(sum);
	a1[5] = u128_high(sum) + a1[5];

	// a1 >>= 1
	a1[0] = (a1[0] >> 1) | (a1[1] << 63);
	a1[1] = (a1[1] >> 1) | (a1[2] << 63);
	a1[2] = (a1[2] >> 1) | (a1[3] << 63);
	a1[3] = (a1[3] >> 1) | (a1[4] << 63);
	a1[4] = (a1[4] >> 1) | (a1[5] << 63);
	a1[5] >>= 1;

	// a1 = (a1 + t) >> (251 = 64 * 3 + 59)
	sum = u128_sum(a1[0], t[0]);
	u128_carry_add(sum, a1[1], t[1]);
	u128_carry_add(sum, a1[2], t[2]);
	u128_carry_add(sum, a1[3], t[3]);
	a1[0] = u128_low(sum);
	u128_carry_add(sum, a1[4], t[4]);
	a1[1] = u128_low(sum);
	a1[2] = u128_high(sum) + a1[5];

	// a1 >>= 59 = a1 / q < 2^(379 - 251 = 128 bits)
	a1[0] = (a1[0] >> 59) | (a1[1] << 5);
	a1[1] = (a1[1] >> 59) | (a1[2] << 5);

	// m' = 0x59D9EBEB3F23782C3164971C4F61FE5CF893F8B602171C89
	static const u64 M2[3] = {
		0xF893F8B602171C89ULL,
		0x3164971C4F61FE5CULL,
		0x59D9EBEB3F23782CULL
	};

	// t <- m' * a2 >> (316 = 64 * 4 + 60)
#define MUL_DIGIT_START(ii, jj) \
	prod = u128_prod_sum(M2[ii], a2[jj], u128_low(sum)); \
	sum = u128_sum(u128_high(sum), u128_high(prod));

#define MUL_DIGIT(ii, jj) \
	prod = u128_prod_sum(M2[ii], a2[jj], u128_low(prod)); \
	u128_add(sum, u128_high(prod));

	prod = u128_prod(M2[0], a2[0]);

	prod = u128_prod_sum(M2[1], a2[0], u128_high(prod));
	u128_set(sum, u128_high(prod));
	MUL_DIGIT(0, 1);

	MUL_DIGIT_START(2, 0);
	MUL_DIGIT(1, 1);
	MUL_DIGIT(0, 2);

	MUL_DIGIT_START(2, 1);
	MUL_DIGIT(1, 2);
	MUL_DIGIT(0, 3);

	MUL_DIGIT_START(2, 2);
	MUL_DIGIT(1, 3);
	MUL_DIGIT(0, 4);

	t[0] = u128_low(prod);

	MUL_DIGIT_START(2, 3);
	MUL_DIGIT(1, 4);

	t[1] = u128_low(prod);

	prod = u128_prod_sum(M2[2], a2[4], u128_low(sum));

	t[2] = u128_low(prod);
	t[3] = u128_high(prod) + u128_high(sum);

#undef MUL_DIGIT_START
#undef MUL_DIGIT

	// t >>= 60
	t[0] = (t[0] >> 60) | (t[1] << 4);
	t[1] = (t[1] >> 60) | (t[2] << 4);
	t[2] = (t[2] >> 60) | (t[3] << 4);
	t[3] >>= 60;

	// a2 -= t
	sum = u128_diff(a2[0], t[0]);
	a2[0] = u128_low(sum);
	u128_borrow_add_sub(sum, a2[1], t[1]);
	a2[1] = u128_low(sum);
	u128_borrow_add_sub(sum, a2[2], t[2]);
	a2[2] = u128_low(sum);
	u128_borrow_add_sub(sum, a2[3], t[3]);
	a2[3] = u128_low(sum);
	a2[4] = u128_high(sum) + a2[4];

	// a2 >>= 1
	a2[0] = (a2[0] >> 1) | (a2[1] << 63);
	a2[1] = (a2[1] >> 1) | (a2[2] << 63);
	a2[2] = (a2[2] >> 1) | (a2[3] << 63);
	a2[3] = (a2[3] >> 1) | (a2[4] << 63);
	a2[4] >>= 1;

	// a2 = (a2 + t) >> (251 = 64 * 3 + 59)
	sum = u128_sum(a2[0], t[0]);
	u128_carry_add(sum, a2[1], t[1]);
	u128_carry_add(sum, a2[2], t[2]);
	u128_carry_add(sum, a2[3], t[3]);
	a2[0] = u128_low(sum);
	a2[1] = u128_high(sum) + a2[4];

	// a2 >>= 59 = a1 / q < 2^(64 bits, verified with MAGMA for worst case k = q - 1)
	a2[0] = (a2[0] >> 59) | (a2[1] << 5);

	// NOTE: The k1, k2 results are guaranteed to fit within 126 bits,
	// so we can safely neglect words 2 and higher since they cannot
	// affect the low 0 and 1 words.

	// t <- (a1 << 126) - a1 = A * a1 < 2^(128+126 = 254)
	// t += B * a2 < 2^(64 + 63 = 127)
	prod = u128_prod(B, a2[0]);
	sum = u128_diff(u128_low(prod), a1[0]);
	t[0] = u128_low(sum);
	t[1] = u128_high(sum) + u128_high(prod) + (a1[0] << 62) - a1[1];

	// k1 <- ||k - t|| = k - (z1 * A + z2 * B) and k1sign
	sum = u128_diff(k[0], t[0]);
	k1.i[0] = u128_low(sum);
	k1.i[1] = u128_high(sum) + k[1] - t[1];
	k1sign = (u32)(k1.i[1] >> 63);

	u64 mask = -(s64)k1sign;

	sum = u128_sum(k1.i[0] ^ mask, k1sign);
	k1.i[0] = u128_low(sum);
	k1.i[1] = u128_high(sum) + (k1.i[1] ^ mask);

	// t <- (a2 << 126) - a2 = A * a2 < 2^(64 + 126 = 190)
	sum = u128_diff((u64)0, a2[0]);
	t[0] = u128_low(sum);
	t[1] = u128_high(sum) + (a2[0] << 62);

	// a1 <- B * a1 < 2^(64+128 = 192)
	prod = u128_prod(B, a1[0]);
	a1[0] = u128_low(prod);
	a1[1] = u128_high(prod) + B * a1[1];

	// k2 <- ||a1 - t|| and k2sign
	sum = u128_diff(a1[0], t[0]);
	k2.i[0] = u128_low(sum);
	u128_borrow_add(sum, a1[1]);
	k2.i[1] = u128_low(sum) - t[1];
	k2sign = (u32)(k2.i[1] >> 63);

	mask = -(s64)k2sign;

	sum = u128_sum(k2.i[0] ^ mask, k2sign);
	k2.i[0] = u128_low(sum);
	k2.i[1] = u128_high(sum) + (k2.i[1] ^ mask);
}

// (x2, y2) = endomorphism(x1, y1)
static CAT_INLINE void gls_morph(const ufe &x1, const ufe &y1, ufe &x2, ufe &y2) {
	// w1 <- conj(x1) = x1^p
	ufe w1;
	fe_conj(x1, w1);

	// x2 <- wx * x1^p
	fe_mul(w1, ENDO_WX, x2);

	// y2 <- y1^p
	fe_conj(y1, y2);
}

// R = endomorphism(P), for extended coordinates
static CAT_INLINE void gls_morph_ext(const ecpt &P, ecpt &R) {
	// w1 <- conj(x1) = x1^p
	ufe w1;
	fe_conj(P.x, w1);

	// x2 <- wx * x1^p
	fe_mul(w1, ENDO_WX, R.x);

	// y2 <- y1^p
	fe_conj(P.y, R.y);

	// For the T coordinate:
	// t = xy/z
	// Z can be neglected.  X and Y are being conjugated, so what is T?
	// t = (a + ib)(c + id)
	// = (ac - bd) + i(bc + ad)
	// t2 = (a - ib)(c - id)
	// = (ac - bd) - i(bc + ad)
	// So t2 = conj(t).

	// w1 <- t1^p
	fe_conj(P.t, w1);

	// x2 <- wx * x1^p
	fe_mul(w1, ENDO_WX, R.t);

	// For the Z coordinate:
	// X1 = x1/z1
	// conj(X1) = conj(x1/z1) = conj(x1) / conj(z1)
	// So the Z coordinate needs to be fixed also

	// z2 <- z1^p
	fe_conj(P.z, R.z);
}

