/*
	Copyright (c) 2013-2014 Christopher A. Taylor.  All rights reserved.

	Redistribution and use in source and binary forms, with or without
	modification, are permitted provided that the following conditions are met:

	* Redistributions of source code must retain the above copyright notice,
	  this list of conditions and the following disclaimer.
	* Redistributions in binary form must reproduce the above copyright notice,
	  this list of conditions and the following disclaimer in the documentation
	  and/or other materials provided with the distribution.
	* Neither the name of Snowshoe nor the names of its contributors may be
	  used to endorse or promote products derived from this software without
	  specific prior written permission.

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
	AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
	ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
	LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
	CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
	SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
	INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
	CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
	ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
	POSSIBILITY OF SUCH DAMAGE.
*/

/*
 * All of the fp_*, fe_*, gls_*, and ec_* math functions:
 * + are branchless and will run in constant time.
 * + are optimized to exploit Fp numbers with high bit = zero.
 * + work regardless of input aliasing (&a == &b == &r is okay).
 * + use incomplete reduction (IR) optimizations where 2^127-1 = 0.
 * + are all in the same large C file for better compiler optimization.
 * + attempt to lower register pressure by using fewer variables.
 * + attempt to take advantage of instruction-level parallelism (ILP).
 * + attempt to make the most of the I-cache with tuned inlining.
 *
 * The input-validation functions are not constant time, but that is not
 * going to leak any information about secret values.  Any functions
 * that are not constant-time are commented with warnings.
 */

#include "ecmul.inc"
#include "snowshoe.h"

#ifndef CAT_ENDIAN_LITTLE

#include "SecureErase.hpp"

/*
 * This file is optimized for little-endian architectures.  In this
 * case the input bytes are already in the internal data format, so
 * it is not necessary to unpack anything.  This is exceptionally
 * good because we need to secure erase any temporary data.  And
 * since no unpacking is done there is no temporary data to erase,
 * and so it will run faster.
 */

static CAT_INLINE void ec_load_k(const char k_chars[32], u64 k[4]) {
	const u64 *k_raw = reinterpret_cast<const u64 *>( k_chars );

	k[0] = getLE(k_raw[0]);
	k[1] = getLE(k_raw[1]);
	k[2] = getLE(k_raw[2]);
	k[3] = getLE(k_raw[3]);
}

static CAT_INLINE void ec_save_k(const u64 k[4], char k_chars[32]) {
	u64 *k_raw = reinterpret_cast<u64 *>( k_chars );

	k_raw[0] = getLE(k[0]);
	k_raw[1] = getLE(k[1]);
	k_raw[2] = getLE(k[2]);
	k_raw[3] = getLE(k[3]);
}

#endif // CAT_ENDIAN_LITTLE

// Check if k == 0 in constant-time
static bool is_zero(const u64 k[4]) {
	u64 zero = k[0] | k[1] | k[2] | k[3];
	u32 z = (u32)(zero | (zero >> 32));
	return z == 0;
}

// Verify that 0 < k < q
static bool invalid_key(const u64 k[4]) {
	// If zero,
	if (is_zero(k)) {
		return true;
	}

	// If not less than q,
	if (!less_q(k)) {
		return true;
	}

	return false;
}

#ifdef __cplusplus
extern "C" {
#endif

int _snowshoe_init(int expected_version) {
	// If math object is aligned oddly.
	if (sizeof(ecpt_affine) != 64) {
		return -1;
	}

	return (expected_version == SNOWSHOE_VERSION) ? 0 : -1;
}

void snowshoe_secret_gen(char k_chars[32]) {
	// Operate on input in-place to avoid making waste variables
	u64 *kq = (u64 *)k_chars;

#ifndef CAT_ENDIAN_LITTLE
	ec_load_k(k_chars, kq);
#endif // CAT_ENDIAN_LITTLE

	ec_mask_scalar(kq);

#ifndef CAT_ENDIAN_LITTLE
	ec_save_k(kq, k_chars);
#endif // CAT_ENDIAN_LITTLE
}

void snowshoe_mul_mod_q(const char x[32], const char y[32], const char z[32], char r[32]) {
#ifndef CAT_ENDIAN_LITTLE
	u64 x1[4+4+4];
	u64 *y1 = x1 + 4;
	u64 *z1 = x1 + 8;

	ec_load_k(x, x1);
	ec_load_k(y, y1);
	if (z) {
		ec_load_k(z, z1);
	}

	mul_mod_q(x1, y1, z ? z1 : 0, x1);

	ec_save_k(x1, r);

	CAT_SECURE_OBJCLR(x1);
#else
	mul_mod_q((const u64 *)x, (const u64 *)y, (const u64 *)z, (u64 *)r);
#endif // CAT_ENDIAN_LITTLE
}

void snowshoe_add_mod_q(const char x[32], const char y[32], char r[32]) {
#ifndef CAT_ENDIAN_LITTLE
	u64 x1[4+4];
	u64 *y1 = x1 + 4;

	ec_load_k(x, x1);
	ec_load_k(y, y1);

	add_mod_q(x1, y1, x1);

	ec_save_k(x1, r);

	CAT_SECURE_OBJCLR(x1);
#else
	add_mod_q((const u64 *)x, (const u64 *)y, (u64 *)r);
#endif // CAT_ENDIAN_LITTLE
}

void snowshoe_mod_q(const char x[64], char r[32]) {
#ifndef CAT_ENDIAN_LITTLE
	u64 x1[8];

	const u64 *k_raw = reinterpret_cast<const u64 *>( x );
	x1[0] = getLE(k_raw[0]);
	x1[1] = getLE(k_raw[1]);
	x1[2] = getLE(k_raw[2]);
	x1[3] = getLE(k_raw[3]);
	x1[4] = getLE(k_raw[4]);
	x1[5] = getLE(k_raw[5]);
	x1[6] = getLE(k_raw[6]);
	x1[7] = getLE(k_raw[7]);

	mod_q(x1, x1);

	ec_save_k(x1, r);

	CAT_SECURE_OBJCLR(x1);
#else
	mod_q((const u64 *)x, (u64 *)r);
#endif // CAT_ENDIAN_LITTLE
}

void snowshoe_neg(const char P[64], char R[64]) {
#ifndef CAT_ENDIAN_LITTLE
	// Load point
	ecpt_affine p1;
	ec_load_xy((const u8 *)P, p1);

	// Run the math routine
	ec_neg_affine(p1, p1);

	// Save result endian-neutral
	ec_save_xy(p1, (u8 *)R);

	CAT_SECURE_OBJCLR(p1); // Maybe unnecessary for all use cases
#else
	ec_neg_affine(*(const ecpt_affine *)P, *(ecpt_affine *)R);
#endif // CAT_ENDIAN_LITTLE
}

int snowshoe_valid(const char P[64]) {
#ifndef CAT_ENDIAN_LITTLE
	// Load point
	ecpt_affine p1;
	ec_load_xy((const u8*)P, p1);

	// If point is invalid,
	if (!ec_valid_vartime(p1)) {
		return -1;
	}

	CAT_SECURE_OBJCLR(p1); // Maybe unnecessary for all use cases
#else
	if (!ec_valid_vartime(*(const ecpt_affine *)P)) {
		return -1;
	}
#endif // CAT_ENDIAN_LITTLE

	return 0;
}

int snowshoe_mul_gen(const char k_raw[32], char R[64], char mul4) {
#ifndef CAT_ENDIAN_LITTLE
	u64 k[4];
	ec_load_k(k_raw, k);

	// Validate key
	if (invalid_key(k)) {
		return -1;
	}

	// R = [4]kG
	ecpt_affine r;
	ecpt p;
	ufe p2b;
	ec_mul_gen(k, p, p2b);
	if (mul4 != 0) {
		ec_dbl(p, p, false, p2b);
		ec_dbl(p, p, false, p2b);
	}
	ec_affine(p, r);

	// Save result endian-neutral
	ec_save_xy(r, (u8*)R);

	CAT_SECURE_OBJCLR(k);
	CAT_SECURE_OBJCLR(r);
#else
	const u64 *k = (const u64 *)k_raw;

	// Validate key
	if (invalid_key(k)) {
		return -1;
	}

	// R = [4]kG
	ecpt p;
	ufe p2b;
	ec_mul_gen(k, p, p2b);
	if (mul4 != 0) {
		ec_dbl(p, p, false, p2b);
		ec_dbl(p, p, false, p2b);
	}
	ec_affine(p, *(ecpt_affine *)R);
#endif // CAT_ENDIAN_LITTLE

	return 0;
}

int snowshoe_mul(const char k_raw[32], const char P[64], char R[64]) {
#ifndef CAT_ENDIAN_LITTLE
	u64 k[4];
	ec_load_k(k_raw, k);

	// Validate key
	if (invalid_key(k)) {
		return -1;
	}

	// Load point
	ecpt_affine p1, r;
	ec_load_xy((const u8*)P, p1);

	// Validate point
	if (!ec_valid(p1)) {
		return -1;
	}

	// Multiply
	ec_mul_affine(k, p1, r);

	// Save result endian-neutral
	ec_save_xy(r, (u8*)R);

	CAT_SECURE_OBJCLR(k);
	CAT_SECURE_OBJCLR(p1);
	CAT_SECURE_OBJCLR(r);
#else
	const u64 *k = (const u64 *)k_raw;

	// Validate key
	if (invalid_key(k)) {
		return -1;
	}

	// Validate point
	if (!ec_valid_vartime(*(const ecpt_affine *)P)) {
		return -1;
	}

	// Multiply
	ec_mul_affine(k, *(const ecpt_affine *)P, *(ecpt_affine *)R);
#endif // CAT_ENDIAN_LITTLE

	return 0;
}

int snowshoe_simul_gen(const char a[32], const char b[32], const char Q[64], char R[64]) {
#ifndef CAT_ENDIAN_LITTLE
	u64 k1[4+4];
	u64 *k2 = k1 + 4;
	ec_load_k(a, k1);
	ec_load_k(b, k2);

	// Validate keys
	if (invalid_key(k1) || invalid_key(k2)) {
		return -1;
	}

	// Load point
	ecpt_affine p2, r;
	ec_load_xy((const u8*)Q, p2);

	// Validate point
	if (!ec_valid(p2)) {
		return -1;
	}

	// Multiply
	ec_simul_gen_affine(k1, k2, p2, r);

	// Save result endian-neutral
	ec_save_xy(r, (u8*)R);

	// May not needed in all use cases
	CAT_SECURE_OBJCLR(k1);
	CAT_SECURE_OBJCLR(p2);
	CAT_SECURE_OBJCLR(r);
#else
	const u64 *k1 = (const u64 *)a;
	const u64 *k2 = (const u64 *)b;
	const ecpt_affine *p2 = (const ecpt_affine *)Q;

	// Validate keys
	if (invalid_key(k1) || invalid_key(k2)) {
		return -1;
	}

	// Validate point
	if (!ec_valid_vartime(*p2)) {
		return -1;
	}

	// Multiply
	ec_simul_gen_affine(k1, k2, *p2, *(ecpt_affine *)R);
#endif // CAT_ENDIAN_LITTLE

	return 0;
}

int snowshoe_simul(const char a[32], const char P[64], const char b[32], const char Q[64], char R[64]) {
#ifndef CAT_ENDIAN_LITTLE
	u64 k1[4], k2[4];
	ec_load_k(a, k1);
	ec_load_k(b, k2);

	// Validate keys
	if (invalid_key(k1) || invalid_key(k2)) {
		return -1;
	}

	// Load points
	ecpt_affine p1, p2, r;
	ec_load_xy((const u8*)P, p1);
	ec_load_xy((const u8*)Q, p2);

	// Validate points
	if (!ec_valid(p1) || !ec_valid(p2)) {
		return -1;
	}

	// Multiply
	ec_simul_affine(k1, p1, k2, p2, r);

	// Save result endian-neutral
	ec_save_xy(r, (u8*)R);

	CAT_SECURE_OBJCLR(k);
	CAT_SECURE_OBJCLR(p1);
	CAT_SECURE_OBJCLR(r);
#else
	const u64 *k1 = (const u64 *)a;
	const u64 *k2 = (const u64 *)b;
	const ecpt_affine *p1 = (const ecpt_affine *)P;
	const ecpt_affine *p2 = (const ecpt_affine *)Q;

	// Validate keys
	if (invalid_key(k1) || invalid_key(k2)) {
		return -1;
	}

	// Validate points
	if (!ec_valid_vartime(*p1) || !ec_valid_vartime(*p2)) {
		return -1;
	}

	// Multiply
	ec_simul_affine(k1, *p1, k2, *p2, *(ecpt_affine *)R);
#endif // CAT_ENDIAN_LITTLE

	return 0;
}

// E = Elligator(key)
int snowshoe_elligator(const char key[32], char E[128]) {
	// Calculate Elligator point from key
	ecpt_affine p;
	ec_elligator_decode(key, p);

	// Validate the resulting point (ie. 0 -> invalid point)
	if (!ec_valid_vartime(p)) {
		return -1;
	}

	// q = 4E
	ecpt q;
	ec_expand(p, q);
	ufe t2b;
	ec_dbl(q, q, true, t2b);
	ec_dbl(q, q, false, t2b);

	// Fix T coordinate
	fe_mul(q.t, t2b, q.t);

	// Copy result
	ecpt *e = (ecpt *)E;
	ec_set(q, *e);

	return 0;
}

// C = kG + E
int snowshoe_elligator_encrypt(const char k[32], const char E[128], char C[64]) {
	// K = kG
	ecpt K;
	const u64 *key = (const u64 *)k;
	if (invalid_key(key)) {
		return -1;
	}
	ufe t2b;
	ec_mul_gen(key, K, t2b);

	// K = K + E
	const ecpt *e = (const ecpt *)E;
	ec_add(K, *e, K, false, false, false, t2b);

	// Affine point
	ecpt_affine *c = (ecpt_affine *)C;
	ec_affine(K, *c);

	return 0;
}

// R = k1(C - E) + k2 * V
int snowshoe_elligator_secret(const char k1[32], const char C[64], const char E[128],
							  const char k2[32], const char V[64], char R[64]) {
	// p = C - E
	ecpt p, q;
	const ecpt_affine *c = (const ecpt_affine *)C;
	if (!ec_valid_vartime(*c)) {
		return -1;
	}
	ec_expand(*c, p);
	const ecpt *e = (const ecpt *)E;
	ec_neg(*e, q);
	ufe t2b;
	ec_add(q, p, p, true, true, true, t2b);

	// If only a single multiplication is required,
	if (!k2) {
		// p = k1 * p
		const u64 *key = (const u64 *)k1;
		if (invalid_key(key)) {
			return -1;
		}
		ec_mul(key, p, false, p, t2b);
	} else {
		// q = V
		const ecpt_affine *v = (const ecpt_affine *)V;
		if (!ec_valid_vartime(*v)) {
			return -1;
		}
		ec_expand(*v, q);

		// p = k1 * p + k2 * q
		const u64 *key1 = (const u64 *)k1;
		const u64 *key2 = (const u64 *)k2;
		ec_simul(key1, p, false, key2, q, true, p, t2b);
	}

	// Fix small subgroup attack
	ec_dbl(p, p, false, t2b);
	ec_dbl(p, p, false, t2b);

	// Affine point
	ecpt_affine *r = (ecpt_affine *)R;
	ec_affine(p, *r);

	return 0;
}

#ifdef __cplusplus
}
#endif

