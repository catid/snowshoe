#include "Platform.hpp"
using namespace cat;

// Prime order of curve
static const u64 EC_Q[4] = {
	0xCE9B68E3B09E01A5ULL,
	0xA6261414C0DC87D3ULL,
	0xFFFFFFFFFFFFFFFFULL,
	0x0FFFFFFFFFFFFFFFULL
};


/*
 * Constant time modular multiplication
 *
d := 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFA6261414C0DC87D3CE9B68E3B09E01A5;

l := 252; // d < 2^l
n := Random(2^256) * Random(2^256);
N := 512;

quot := n div d;
rem := n mod d;

mp := 2^(N+l) div d - 2^N + 1;

print mp;
print 221269984318908782116022557274131459378616010318126502992777998111130512380382470982435267452138673041816973846238807:Hex;

// Live computation starts here:

// t <- MULUH m' * n
t := (mp * n) div (2^N);

// s <- t + SRL(n - t, 1)
s := t + ((n - t) div 2);

// quotient <- s >> (l - 1)
quotp := s div (2^(l-1));

print quot;
print quotp;

remp := n - quot * d;

print remp;
print rem;
 * 
 */

// r = x * y + z (mod q), z optional
static void mul_mod_q(const u64 x[4], const u64 y[4], const u64 z[4], u64 r[4]) {
	u64 p[8], t[7], n[4];
	s128 diff;
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

	// n = p (low part), needed at end
	n[0] = p[0];
	n[1] = p[1];
	n[2] = p[2];
	n[3] = p[3];

	/*
	 * Using the Unsigned Division algorithm from section 4 of [2]:
	 *
	 * Computing quot = n / d
	 *
	 * d = q < 2^252, l = 252
	 * p < 2^(252+256), N = 252+256
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
	 * See magma_unsigned_division.txt for more details.
	 */

	// m' = 0x59D9EBEB3F23782C3164971C4F61FE5CF893F8B602171C88E95EB7B0E1A988566D91A79575334CACB91DD2622FBD3D66
	static const u64 M1[6] = {
		0xB91DD2622FBD3D66ULL,
		0x6D91A79575334CACULL,
		0xE95EB7B0E1A98856ULL,
		0xF893F8B602171C88ULL,
		0x3164971C4F61FE5CULL,
		0x59D9EBEB3F23782CULL
	};

	// t <- m' * a1 >> (252 + 256 = 508 = 64*7 + 60)

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

	prod = (u128)M1[5] * p[1] + (u64)carry;
	carry >>= 64;
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

	prod = (u128)M1[5] * p[2] + (u64)carry;
	carry >>= 64;
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

	t[0] = (u64)prod;

	prod = (u128)M1[5] * p[3] + (u64)carry;
	carry >>= 64;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[4] * p[4] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[3] * p[5] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[2] * p[6] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[1] * p[7] + (u64)prod;
	carry += (u64)(prod >> 64);

	t[1] = (u64)prod;

	prod = (u128)M1[5] * p[4] + (u64)carry;
	carry >>= 64;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[4] * p[5] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[3] * p[6] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[2] * p[7] + (u64)prod;
	carry += (u64)(prod >> 64);

	t[2] = (u64)prod;

	prod = (u128)M1[5] * p[5] + (u64)carry;
	carry >>= 64;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[4] * p[6] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[3] * p[7] + (u64)prod;
	carry += (u64)(prod >> 64);

	t[3] = (u64)prod;

	prod = (u128)M1[5] * p[6] + (u64)carry;
	carry >>= 64;
	carry += (u64)(prod >> 64);
	prod = (u128)M1[4] * p[7] + (u64)prod;
	carry += (u64)(prod >> 64);

	t[4] = (u64)prod;

	prod = (u128)M1[5] * p[7] + (u64)carry;
	carry >>= 64;
	carry += (u64)(prod >> 64);

	t[5] = (u64)prod;
	t[6] = (u64)carry;

	// t >>= 60
	t[0] = (t[0] >> 60) | (t[1] << 4);
	t[1] = (t[1] >> 60) | (t[2] << 4);
	t[2] = (t[2] >> 60) | (t[3] << 4);
	t[3] = (t[3] >> 60) | (t[4] << 4);
	t[4] = (t[4] >> 60) | (t[5] << 4);
	t[5] = (t[5] >> 60) | (t[6] << 4);
	t[6] >>= 60;

	// p -= t
	diff = p[0] - t[0];
	p[0] = (u64)diff;
	diff = ((diff >> 64) + p[1]) - t[1];
	p[1] = (u64)diff;
	diff = ((diff >> 64) + p[2]) - t[2];
	p[2] = (u64)diff;
	diff = ((diff >> 64) + p[3]) - t[3];
	p[3] = (u64)diff;
	diff = ((diff >> 64) + p[4]) - t[4];
	p[4] = (u64)diff;
	diff = ((diff >> 64) + p[5]) - t[5];
	p[5] = (u64)diff;
	diff = ((diff >> 64) + p[6]) - t[6];
	p[6] = (u64)diff;
	diff = (diff >> 64) + p[7];
	p[7] = (u64)diff;

	// p >>= 1
	p[0] = (p[0] >> 1) | (p[1] << 63);
	p[1] = (p[1] >> 1) | (p[2] << 63);
	p[2] = (p[2] >> 1) | (p[3] << 63);
	p[3] = (p[3] >> 1) | (p[4] << 63);
	p[4] = (p[4] >> 1) | (p[5] << 63);
	p[5] = (p[5] >> 1) | (p[6] << 63);
	p[6] = (p[6] >> 1) | (p[7] << 63);
	p[7] >>= 1;

	// p = (p + t) >> (251 = 64 * 3 + 59)
	sum = (u128)p[0] + t[0];
	sum = ((u128)p[1] + t[1]) + (u64)(sum >> 64);
	sum = ((u128)p[2] + t[2]) + (u64)(sum >> 64);
	sum = ((u128)p[3] + t[3]) + (u64)(sum >> 64);
	p[0] = (u64)sum;
	sum = ((u128)p[4] + t[4]) + (u64)(sum >> 64);
	p[1] = (u64)sum;
	sum = ((u128)p[5] + t[5]) + (u64)(sum >> 64);
	p[2] = (u64)sum;
	sum = ((u128)p[6] + t[6]) + (u64)(sum >> 64);
	p[3] = (u64)sum;
	sum = (u128)p[7] + (u64)(sum >> 64);
	p[4] = (u64)sum;

	// p >>= 59 = a1 / q < 2^(252 + 256 - 251 = 257 bits)
	p[0] = (p[0] >> 59) | (p[1] << 5);
	p[1] = (p[1] >> 59) | (p[2] << 5);
	p[2] = (p[2] >> 59) | (p[3] << 5);
	p[3] = (p[3] >> 59) | (p[4] << 5);

	// NOTE: p is now the quotient of (x * y + z) / q
	// To recover the remainder, we need to multiply by q again:

	// p = p * q (only need low 4 words of it)

	// Comba multiplication: Right to left schoolbook column approach
	prod = (u128)p[0] * EC_Q[0];
	p[0] = (u64)prod;

	prod = (u128)p[1] * EC_Q[0] + (u64)(prod >> 64);
	carry = (u64)(prod >> 64);
	prod = (u128)p[0] * EC_Q[1] + (u64)prod;
	carry += (u64)(prod >> 64);
	p[1] = (u64)prod;

	prod = (u128)p[2] * EC_Q[0] + (u64)carry;
	carry >>= 64;
	carry += (u64)(prod >> 64);
	prod = (u128)p[1] * EC_Q[1] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)p[0] * EC_Q[2] + (u64)prod;
	carry += (u64)(prod >> 64);
	p[2] = (u64)prod;

	prod = (u128)p[3] * EC_Q[0] + (u64)carry;
	carry >>= 64;
	carry += (u64)(prod >> 64);
	prod = (u128)p[2] * EC_Q[1] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)p[1] * EC_Q[2] + (u64)prod;
	carry += (u64)(prod >> 64);
	prod = (u128)p[0] * EC_Q[3] + (u64)prod;
	p[3] = (u64)prod;

	// And then subtract it from the original result to get the remainder:

	// r = n - p
	diff = n[0] - p[0];
	r[0] = (u64)diff;
	diff = ((diff >> 64) + n[1]) - p[1];
	r[1] = (u64)diff;
	diff = ((diff >> 64) + n[2]) - p[2];
	r[2] = (u64)diff;
	diff = ((diff >> 64) + n[3]) - p[3];
	r[3] = (u64)diff;
}

