// Prime order of curve
static const u64 EC_Q[4] = {
	0xCE9B68E3B09E01A5ULL,
	0xA6261414C0DC87D3ULL,
	0xFFFFFFFFFFFFFFFFULL,
	0x0FFFFFFFFFFFFFFFULL
};

/*
 * Mask a random number to produce a compatible scalar for multiplication
 */

static CAT_INLINE void ec_mask_scalar(u64 k[4]) {
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

// Returns true if k is less than q
static CAT_INLINE bool less_q(const u64 k[4]) {
	u128 sum;

	// p -= t
	sum = u128_diff(k[0], EC_Q[0]);
	u128_borrow_add_sub(sum, k[1], EC_Q[1]);
	u128_borrow_add_sub(sum, k[2], EC_Q[2]);
	u128_borrow_add_sub(sum, k[3], EC_Q[3]);

	// If there was a borrow out, then it was less than q
	return u128_is_neg(sum);
}

// r = q - k
static CAT_INLINE void neg_mod_q(const u64 k[4], u64 r[4]) {
	u128 sum;

	// p -= t
	sum = u128_diff(EC_Q[0], k[0]);
	r[0] = u128_low(sum);
	u128_borrow_add_sub(sum, EC_Q[1], k[1]);
	r[1] = u128_low(sum);
	u128_borrow_add_sub(sum, EC_Q[2], k[2]);
	r[2] = u128_low(sum);
	u128_borrow_add_sub(sum, EC_Q[3], k[3]);
	r[3] = u128_low(sum);
}

// r = p (mod q), where p < 2^512
static void mod_q(const u64 p[8], u64 r[4]) {
	u64 t[7], n[8];
	u128 sum, prod;

	/*
	 * Using the Unsigned Division algorithm from section 4 of [2]:
	 *
	 * Computing quot = n / d
	 *
	 * d = q < 2^252, l = 252
	 * p < 2^512, N = 512
	 *
	 * m' = 2^N * (2^l - d) / d + 1
	 * = 2^(N+l)/d - 2^N + 1
	 *
	 * t := (m' * p) div (2^N);
	 *
	 * s := t + ((p - t) div 2);
	 *
	 * quot := s div (2^(l-1));
	 *
	 * See magma_unsigned_remainder.txt for more details.
	 */

	// m' = 0x59D9EBEB3F23782C3164971C4F61FE5CF893F8B602171C88E95EB7B0E1A988566D91A79575334CACB91DD2622FBD3D657
	static const u64 M1[7] = {
		0x91DD2622FBD3D657ULL,
		0xD91A79575334CACBULL,
		0x95EB7B0E1A988566ULL,
		0x893F8B602171C88EULL,
		0x164971C4F61FE5CFULL,
		0x9D9EBEB3F23782C3ULL,
		5 // replace multiplications by M1[6] with shift+add
	};

	// t <- m' * p >> (512 = 64*8)

	// Comba multiplication: Right to left schoolbook column approach
#define MUL_DIGIT_START(ii, jj) \
	prod = u128_prod_sum(M1[ii], p[jj], u128_low(sum)); \
	sum = u128_sum(u128_high(sum), u128_high(prod));

#define MUL_DIGIT(ii, jj) \
	prod = u128_prod_sum(M1[ii], p[jj], u128_low(prod)); \
	u128_add(sum, u128_high(prod));

	prod = u128_prod(M1[0], p[0]);

	prod = u128_prod_sum(M1[1], p[0], u128_high(prod));
	u128_set(sum, u128_high(prod));
	MUL_DIGIT(0, 1);

	MUL_DIGIT_START(2, 0);
	MUL_DIGIT(1, 1);
	MUL_DIGIT(0, 2);

	MUL_DIGIT_START(3, 0);
	MUL_DIGIT(2, 1);
	MUL_DIGIT(1, 2);
	MUL_DIGIT(0, 3);

	MUL_DIGIT_START(4, 0);
	MUL_DIGIT(3, 1);
	MUL_DIGIT(2, 2);
	MUL_DIGIT(1, 3);
	MUL_DIGIT(0, 4);

	MUL_DIGIT_START(5, 0);
	MUL_DIGIT(4, 1);
	MUL_DIGIT(3, 2);
	MUL_DIGIT(2, 3);
	MUL_DIGIT(1, 4);
	MUL_DIGIT(0, 5);

#undef MUL_DIGIT_START

#define MUL_DIGIT_START(jj) \
	prod = u128_lshift_sum(p[jj], 2, p[jj]); \
	u128_add(prod, u128_low(sum)); \
	sum = u128_sum(u128_high(sum), u128_high(prod));

	MUL_DIGIT_START(0);
	MUL_DIGIT(5, 1);
	MUL_DIGIT(4, 2);
	MUL_DIGIT(3, 3);
	MUL_DIGIT(2, 4);
	MUL_DIGIT(1, 5);
	MUL_DIGIT(0, 6);

	MUL_DIGIT_START(1);
	MUL_DIGIT(5, 2);
	MUL_DIGIT(4, 3);
	MUL_DIGIT(3, 4);
	MUL_DIGIT(2, 5);
	MUL_DIGIT(1, 6);
	MUL_DIGIT(0, 7);

	MUL_DIGIT_START(2);
	MUL_DIGIT(5, 3);
	MUL_DIGIT(4, 4);
	MUL_DIGIT(3, 5);
	MUL_DIGIT(2, 6);
	MUL_DIGIT(1, 7);

	t[0] = u128_low(prod);

	MUL_DIGIT_START(3);
	MUL_DIGIT(5, 4);
	MUL_DIGIT(4, 5);
	MUL_DIGIT(3, 6);
	MUL_DIGIT(2, 7);

	t[1] = u128_low(prod);

	MUL_DIGIT_START(4);
	MUL_DIGIT(5, 5);
	MUL_DIGIT(4, 6);
	MUL_DIGIT(3, 7);

	t[2] = u128_low(prod);

	MUL_DIGIT_START(5);
	MUL_DIGIT(5, 6);
	MUL_DIGIT(4, 7);

	t[3] = u128_low(prod);

	MUL_DIGIT_START(6);
	MUL_DIGIT(5, 7);

	t[4] = u128_low(prod);

	prod = u128_lshift_sum(p[7], 2, p[7]);
	u128_add(prod, u128_low(sum));

	t[5] = u128_low(prod);
	t[6] = u128_high(prod) + u128_high(sum);

#undef MUL_DIGIT_START
#undef MUL_DIGIT

	// n = p - t
	sum = u128_diff(p[0], t[0]);
	n[0] = u128_low(sum);
	u128_borrow_add_sub(sum, p[1], t[1]);
	n[1] = u128_low(sum);
	u128_borrow_add_sub(sum, p[2], t[2]);
	n[2] = u128_low(sum);
	u128_borrow_add_sub(sum, p[3], t[3]);
	n[3] = u128_low(sum);
	u128_borrow_add_sub(sum, p[4], t[4]);
	n[4] = u128_low(sum);
	u128_borrow_add_sub(sum, p[5], t[5]);
	n[5] = u128_low(sum);
	u128_borrow_add_sub(sum, p[6], t[6]);
	n[6] = u128_low(sum);
	n[7] = u128_high(sum) + p[7];

	// n >>= 1
	n[0] = (n[0] >> 1) | (n[1] << 63);
	n[1] = (n[1] >> 1) | (n[2] << 63);
	n[2] = (n[2] >> 1) | (n[3] << 63);
	n[3] = (n[3] >> 1) | (n[4] << 63);
	n[4] = (n[4] >> 1) | (n[5] << 63);
	n[5] = (n[5] >> 1) | (n[6] << 63);
	n[6] = (n[6] >> 1) | (n[7] << 63);
	n[7] >>= 1;

	// n = (n + t) >> (251 = 64 * 3 + 59)
	sum = u128_sum(n[0], t[0]);
	u128_carry_add(sum, n[1], t[1]);
	u128_carry_add(sum, n[2], t[2]);
	u128_carry_add(sum, n[3], t[3]);
	n[0] = u128_low(sum);
	u128_carry_add(sum, n[4], t[4]);
	n[1] = u128_low(sum);
	u128_carry_add(sum, n[5], t[5]);
	n[2] = u128_low(sum);
	u128_carry_add(sum, n[6], t[6]);
	n[3] = u128_low(sum);
	n[4] = u128_high(sum) + n[7];

	// n >>= 59 = p / q < 2^(512 - 251 = 261 bits)
	n[0] = (n[0] >> 59) | (n[1] << 5);
	n[1] = (n[1] >> 59) | (n[2] << 5);
	n[2] = (n[2] >> 59) | (n[3] << 5);
	n[3] = (n[3] >> 59) | (n[4] << 5);
	//n[4] = n[4] >> 59;

	// NOTE: n is now the quotient of p / q
	// To recover the remainder, we need to multiply by q again:

	// t = n * q (only need low 4 words of it)

	// Comba multiplication: Right to left schoolbook column approach
	prod = u128_prod(n[0], EC_Q[0]);

	t[0] = u128_low(prod);

	prod = u128_prod_sum(n[1], EC_Q[0], u128_high(prod));
	u128_set(sum, u128_high(prod));
	prod = u128_prod_sum(n[0], EC_Q[1], u128_low(prod));
	u128_add(sum, u128_high(prod));

	t[1] = u128_low(prod);

	prod = u128_prod_sum(n[2], EC_Q[0], u128_low(sum));
	u64 temp = u128_high(sum) + u128_high(prod);
	prod = u128_prod_sum(n[1], EC_Q[1], u128_low(prod));
	temp += u128_high(prod);
	prod = u128_prod_sum(n[0], EC_Q[2], u128_low(prod));
	temp += u128_high(prod);

	t[2] = u128_low(prod);
	t[3] = temp + n[3] * EC_Q[0] + n[2] * EC_Q[1] + n[1] * EC_Q[2] + n[0] * EC_Q[3];

	// And then subtract it from the original input to get the remainder:

	// r = p - t
	sum = u128_diff(p[0], t[0]);
	r[0] = u128_low(sum);
	u128_borrow_add_sub(sum, p[1], t[1]);
	r[1] = u128_low(sum);
	u128_borrow_add_sub(sum, p[2], t[2]);
	r[2] = u128_low(sum);
	r[3] = u128_high(sum) + p[3] - t[3];
}

// r = x * y + z (mod q), z optional
static void mul_mod_q(const u64 x[4], const u64 y[4], const u64 z[4], u64 r[4]) {
	u64 p[8];
	u128 sum, prod;

	// p = x * y < 2^(256 + 252 = 508)

	// Comba multiplication: Right to left schoolbook column approach
#define MUL_DIGIT_START(ii, jj) \
	prod = u128_prod_sum(x[ii], y[jj], u128_low(sum)); \
	sum = u128_sum(u128_high(sum), u128_high(prod));

#define MUL_DIGIT(ii, jj) \
	prod = u128_prod_sum(x[ii], y[jj], u128_low(prod)); \
	u128_add(sum, u128_high(prod));

	prod = u128_prod(x[0], y[0]);

	p[0] = u128_low(prod);

	prod = u128_prod_sum(x[1], y[0], u128_high(prod));
	u128_set(sum, u128_high(prod));
	MUL_DIGIT(0, 1);

	p[1] = u128_low(prod);

	MUL_DIGIT_START(2, 0);
	MUL_DIGIT(1, 1);
	MUL_DIGIT(0, 2);

	p[2] = u128_low(prod);

	MUL_DIGIT_START(3, 0);
	MUL_DIGIT(2, 1);
	MUL_DIGIT(1, 2);
	MUL_DIGIT(0, 3);

	p[3] = u128_low(prod);

	MUL_DIGIT_START(3, 1);
	MUL_DIGIT(2, 2);
	MUL_DIGIT(1, 3);

	p[4] = u128_low(prod);

	MUL_DIGIT_START(3, 2);
	MUL_DIGIT(2, 3);

	p[5] = u128_low(prod);

	prod = u128_prod_sum(x[3], y[3], u128_low(sum));

	p[6] = u128_low(prod);
	p[7] = u128_high(prod) + u128_high(sum);

#undef MUL_DIGIT_START
#undef MUL_DIGIT

	// If z is provided,
	if (z) {
		sum = u128_sum(p[0], z[0]);
		p[0] = u128_low(sum);
		u128_carry_add(sum, p[1], z[1]);
		p[1] = u128_low(sum);
		u128_carry_add(sum, p[2], z[2]);
		p[2] = u128_low(sum);
		u128_carry_add(sum, p[3], z[3]);
		p[3] = u128_low(sum);
		u128_carry_add(sum, p[4]);
		p[4] = u128_low(sum);
		u128_carry_add(sum, p[5]);
		p[5] = u128_low(sum);
		u128_carry_add(sum, p[6]);
		p[6] = u128_low(sum);
		p[7] = u128_high(sum) + p[7];
	}

	mod_q(p, r);
}

// r = x + y (mod q)
static CAT_INLINE void add_mod_q(const u64 x[4], const u64 y[4], u64 r[4]) {
	u64 z[8];
	u128 sum;

	// z = x + y
	sum = u128_sum(x[0], y[0]);
	z[0] = u128_low(sum);
	u128_carry_add(sum, x[1], y[1]);
	z[1] = u128_low(sum);
	u128_carry_add(sum, x[2], y[2]);
	z[2] = u128_low(sum);
	u128_carry_add(sum, x[3], y[3]);
	z[3] = u128_low(sum);

	// Pad the result out to 512 bits, so that it can be used by mod_q()
	z[4] = u128_high(sum);
	z[5] = 0;
	z[6] = 0;
	z[7] = 0;

	// Note this is not incredibly efficient, but I would rather only have to
	// audit one version of mod_q instead of two.

	mod_q(z, r);
}

