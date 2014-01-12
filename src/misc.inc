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
	// p -= t
	s128 diff = (s128)k[0] - EC_Q[0];
	diff = ((diff >> 64) + k[1]) - EC_Q[1];
	diff = ((diff >> 64) + k[2]) - EC_Q[2];
	diff = ((diff >> 64) + k[3]) - EC_Q[3];

	// If there was a borrow out, then it was less than q
	return diff < 0;
}

// r = q - k
static CAT_INLINE void neg_mod_q(const u64 k[4], u64 r[4]) {
	s128 diff = (s128)EC_Q[0] - k[0];
	r[0] = (u64)diff;
	diff = ((diff >> 64) + EC_Q[1]) - k[1];
	r[1] = (u64)diff;
	diff = ((diff >> 64) + EC_Q[2]) - k[2];
	r[2] = (u64)diff;
	diff = ((diff >> 64) + EC_Q[3]) - k[3];
	r[3] = (u64)diff;
}

// r = p (mod q), where p < 2^512
static void mod_q(const u64 p[8], u64 r[4]) {
	u64 t[7], n[8];
	s128 diff;
	u128 prod, sum, carry;

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
	prod = (u128)M1[0] * p[0];

	prod = (u128)M1[1] * p[0] + (u64)(prod >> 64);
	carry = (u64)(prod >> 64);
	prod = (u128)M1[0] * p[1] + (u64)prod;
	carry += (u64)(prod >> 64);

	prod = (u128)M1[2] * p[0] + (u64)carry;
	carry >>= 64;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[1] * p[1] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[0] * p[2] + (u64)prod;
	carry += (u64)(prod >> 64);

	prod = (u128)M1[3] * p[0] + (u64)carry;
	carry >>= 64;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[2] * p[1] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[1] * p[2] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[0] * p[3] + (u64)prod;
	carry += (u64)(prod >> 64);

	prod = (u128)M1[4] * p[0] + (u64)carry;
	carry >>= 64;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[3] * p[1] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[2] * p[2] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[1] * p[3] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[0] * p[4] + (u64)prod;
	carry += (u64)(prod >> 64);

	prod = (u128)M1[5] * p[0] + (u64)carry;
	carry >>= 64;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[4] * p[1] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[3] * p[2] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[2] * p[3] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[1] * p[4] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[0] * p[5] + (u64)prod;
	carry += (u64)(prod >> 64);

	prod = ((u128)p[0] << 2) + p[0] + (u64)carry;
	carry >>= 64;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[5] * p[1] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[4] * p[2] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[3] * p[3] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[2] * p[4] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[1] * p[5] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[0] * p[6] + (u64)prod;
	carry += (u64)(prod >> 64);

	prod = ((u128)p[1] << 2) + p[1] + (u64)carry;
	carry >>= 64;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[5] * p[2] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[4] * p[3] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[3] * p[4] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[2] * p[5] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[1] * p[6] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[0] * p[7] + (u64)prod;
	carry += (u64)(prod >> 64);

	prod = ((u128)p[2] << 2) + p[2] + (u64)carry;
	carry >>= 64;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[5] * p[3] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[4] * p[4] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[3] * p[5] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[2] * p[6] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[1] * p[7] + (u64)prod;
	carry += (u64)(prod >> 64);

	t[0] = (u64)prod;

	prod = ((u128)p[3] << 2) + p[3] + (u64)carry;
	carry >>= 64;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[5] * p[4] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[4] * p[5] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[3] * p[6] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[2] * p[7] + (u64)prod;
	carry += (u64)(prod >> 64);

	t[1] = (u64)prod;

	prod = ((u128)p[4] << 2) + p[4] + (u64)carry;
	carry >>= 64;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[5] * p[5] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[4] * p[6] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[3] * p[7] + (u64)prod;
	carry += (u64)(prod >> 64);

	t[2] = (u64)prod;

	prod = ((u128)p[5] << 2) + p[5] + (u64)carry;
	carry >>= 64;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[5] * p[6] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[4] * p[7] + (u64)prod;
	carry += (u64)(prod >> 64);

	t[3] = (u64)prod;

	prod = ((u128)p[6] << 2) + p[6] + (u64)carry;
	carry >>= 64;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[5] * p[7] + (u64)prod;
	carry += (u64)(prod >> 64);

	t[4] = (u64)prod;

	prod = ((u128)p[7] << 2) + p[7] + (u64)carry;
	carry >>= 64;
	carry += (u64)(prod >> 64);

	t[5] = (u64)prod;
	t[6] = (u64)carry;

	// n = p - t
	diff = (s128)p[0] - t[0];
	n[0] = (u64)diff;
	diff = ((diff >> 64) + p[1]) - t[1];
	n[1] = (u64)diff;
	diff = ((diff >> 64) + p[2]) - t[2];
	n[2] = (u64)diff;
	diff = ((diff >> 64) + p[3]) - t[3];
	n[3] = (u64)diff;
	diff = ((diff >> 64) + p[4]) - t[4];
	n[4] = (u64)diff;
	diff = ((diff >> 64) + p[5]) - t[5];
	n[5] = (u64)diff;
	diff = ((diff >> 64) + p[6]) - t[6];
	n[6] = (u64)diff;
	diff = (diff >> 64) + p[7];
	n[7] = (u64)diff;

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
	sum = (u128)n[0] + t[0];
	sum = ((u128)n[1] + t[1]) + (u64)(sum >> 64);
	sum = ((u128)n[2] + t[2]) + (u64)(sum >> 64);
	sum = ((u128)n[3] + t[3]) + (u64)(sum >> 64);
	n[0] = (u64)sum;
	sum = ((u128)n[4] + t[4]) + (u64)(sum >> 64);
	n[1] = (u64)sum;
	sum = ((u128)n[5] + t[5]) + (u64)(sum >> 64);
	n[2] = (u64)sum;
	sum = ((u128)n[6] + t[6]) + (u64)(sum >> 64);
	n[3] = (u64)sum;
	sum = (u128)n[7] + (u64)(sum >> 64);
	n[4] = (u64)sum;

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
	prod = (u128)n[0] * EC_Q[0];
	t[0] = (u64)prod;

	prod = (u128)n[1] * EC_Q[0] + (u64)(prod >> 64);
	carry = (u64)(prod >> 64);
	prod = (u128)n[0] * EC_Q[1] + (u64)prod;
	carry += (u64)(prod >> 64);
	t[1] = (u64)prod;

	prod = (u128)n[2] * EC_Q[0] + (u64)carry;
	carry >>= 64;
	carry += (u64)(prod >> 64);
	prod = (u128)n[1] * EC_Q[1] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)n[0] * EC_Q[2] + (u64)prod;
	carry += (u64)(prod >> 64);
	t[2] = (u64)prod;

	prod = (u128)n[3] * EC_Q[0] + (u64)carry;
	carry >>= 64;
	carry += (u64)(prod >> 64);
	prod = (u128)n[2] * EC_Q[1] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)n[1] * EC_Q[2] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)n[0] * EC_Q[3] + (u64)prod;
	t[3] = (u64)prod;

	// And then subtract it from the original input to get the remainder:

	// r = p - t
	diff = (s128)p[0] - t[0];
	r[0] = (u64)diff;
	diff = ((diff >> 64) + p[1]) - t[1];
	r[1] = (u64)diff;
	diff = ((diff >> 64) + p[2]) - t[2];
	r[2] = (u64)diff;
	diff = ((diff >> 64) + p[3]) - t[3];
	r[3] = (u64)diff;
}

// r = x * y + z (mod q), z optional
static void mul_mod_q(const u64 x[4], const u64 y[4], const u64 z[4], u64 r[4]) {
	u64 p[8];
	u128 prod, sum, carry;

	// p = x * y < 2^(256 + 252 = 508)

	// Comba multiplication: Right to left schoolbook column approach
	prod = (u128)x[0] * y[0];
	p[0] = (u64)prod;

	prod = (u128)x[1] * y[0] + (u64)(prod >> 64);
	carry = (u64)(prod >> 64);
	prod = (u128)x[0] * y[1] + (u64)prod;
	carry += (u64)(prod >> 64);
	p[1] = (u64)prod;

	prod = (u128)x[2] * y[0] + (u64)carry;
	carry >>= 64;
	carry += (u64)(prod >> 64);
	prod = (u128)x[1] * y[1] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)x[0] * y[2] + (u64)prod;
	carry += (u64)(prod >> 64);
	p[2] = (u64)prod;

	prod = (u128)x[3] * y[0] + (u64)carry;
	carry >>= 64;
	carry += (u64)(prod >> 64);
	prod = (u128)x[2] * y[1] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)x[1] * y[2] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)x[0] * y[3] + (u64)prod;
	carry += (u64)(prod >> 64);
	p[3] = (u64)prod;

	prod = (u128)x[3] * y[1] + (u64)carry;
	carry >>= 64;
	carry += (u64)(prod >> 64);
	prod = (u128)x[2] * y[2] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)x[1] * y[3] + (u64)prod;
	carry += (u64)(prod >> 64);
	p[4] = (u64)prod;

	prod = (u128)x[3] * y[2] + (u64)carry;
	carry >>= 64;
	carry += (u64)(prod >> 64);
	prod = (u128)x[2] * y[3] + (u64)prod;
	carry += (u64)(prod >> 64);
	p[5] = (u64)prod;

	prod = (u128)x[3] * y[3] + (u64)carry;
	carry >>= 64;
	carry += (u64)(prod >> 64);
	p[6] = (u64)prod;
	p[7] = (u64)carry;

	// If z is provided,
	if (z) {
		sum = (u128)p[0] + z[0];
		p[0] = (u64)sum;
		sum = ((u128)p[1] + z[1]) + (u64)(sum >> 64);
		p[1] = (u64)sum;
		sum = ((u128)p[2] + z[2]) + (u64)(sum >> 64);
		p[2] = (u64)sum;
		sum = ((u128)p[3] + z[3]) + (u64)(sum >> 64);
		p[3] = (u64)sum;
		sum = (u128)p[4] + (u64)(sum >> 64);
		p[4] = (u64)sum;
		sum = (u128)p[5] + (u64)(sum >> 64);
		p[5] = (u64)sum;
		sum = (u128)p[6] + (u64)(sum >> 64);
		p[6] = (u64)sum;
		sum = (u128)p[7] + (u64)(sum >> 64);
		p[7] = (u64)sum;
	}

	mod_q(p, r);
}

// r = x + y (mod q)
static CAT_INLINE void add_mod_q(const u64 x[4], const u64 y[4], u64 r[4]) {
	u64 z[8];

	// z = x + y
	u128 sum = (u128)x[0] + y[0];
	z[0] = (u64)sum;
	sum = ((u128)x[1] + y[1]) + (u64)(sum >> 64);
	z[1] = (u64)sum;
	sum = ((u128)x[2] + y[2]) + (u64)(sum >> 64);
	z[2] = (u64)sum;
	sum = ((u128)x[3] + y[3]) + (u64)(sum >> 64);
	z[3] = (u64)sum;

	// Pad the result out to 512 bits, so that it can be used by mod_q()
	z[4] = (u64)(sum >> 64);
	z[5] = 0;
	z[6] = 0;
	z[7] = 0;

	// Note this is not incredibly efficient, but I would rather only have to
	// audit one version of mod_q instead of two.

	mod_q(z, r);
}

