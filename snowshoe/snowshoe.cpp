/*
	Copyright (c) 2013 Christopher A. Taylor.  All rights reserved.

	Redistribution and use in source and binary forms, with or without
	modification, are permitted provided that the following conditions are met:

	* Redistributions of source code must retain the above copyright notice,
	  this list of conditions and the following disclaimer.
	* Redistributions in binary form must reproduce the above copyright notice,
	  this list of conditions and the following disclaimer in the documentation
	  and/or other materials provided with the distribution.
	* Neither the name of Tabby nor the names of its contributors may be
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
 * + expect the input to have a clear high bit.
 * + work regardless of input aliasing (&a == &b == &r is okay).
 * + are inline to allow for full compiler optimization.
 * + attempt to use as few registers as possible.
 * + attempt to take advantage of CPU pipelining.
 *
 * The input-validation functions are not constant time, but that is not
 * going to leak any important information about long-term secrets.  Any
 * functions that are not constant-time are commented with warnings.
 */

#include "ecmul.inc"
#include "snowshoe.h"

#ifdef __cplusplus
extern "C" {
#endif

bool _snowshoe_init(int expected_version) {
	return expected_version == SNOWSHOE_VERSION;
}

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

void snowshoe_secret_gen(char k_chars[32]) {
	u64 kq[4];
	ec_load_k(k_chars, kq);

	ec_mask_scalar(kq);

	ec_save_k(kq, k_chars);
}

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

void snowshoe_mul_mod_q(const char x[32], const char y[32], const char z[32], char r[32]) {
	u64 x1[4], y1[4], z1[4];
	ec_load_k(x, x1);
	ec_load_k(y, y1);
	if (z) {
		ec_load_k(z, z1);
	}

	mul_mod_q(x1, y1, z ? z1 : 0, x1);

	ec_save_k(x1, r);
}

void snowshoe_mod_q(const char x[64], char r[32]) {
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
}

void snowshoe_neg(const char P[64], char R[64]) {
	// Load point
	ecpt_affine p1;
	ec_load_xy((const u8*)P, p1);

	// Run the math routine
	ec_neg_affine(p1, p1);

	// Save result endian-neutral
	ec_save_xy(p1, (u8*)R);
}

bool snowshoe_mul_gen(const char k_raw[32], const bool mul_cofactor, const bool constant_time, char R[64]) {
	u64 k[4];
	ec_load_k(k_raw, k);

	// Validate key
	if (invalid_key(k)) {
		return false;
	}

	// Run the math routine
	ecpt_affine r;
	ec_mul_gen(k, mul_cofactor, constant_time, r);

	// Save result endian-neutral
	ec_save_xy(r, (u8*)R);

	return true;
}

bool snowshoe_mul(const char k_raw[32], const char P[64], char R[64]) {
	u64 k[4];
	ec_load_k(k_raw, k);

	// Validate key
	if (invalid_key(k)) {
		return false;
	}

	// Load point
	ecpt_affine p1, r;
	ec_load_xy((const u8*)P, p1);

	// Validate point
	if (!ec_valid(p1.x, p1.y)) {
		return false;
	}

	// Run the math routine
	ec_mul(k, p1, r);

	// Save result endian-neutral
	ec_save_xy(r, (u8*)R);

	return true;
}

bool snowshoe_simul_gen(const char a[32], const char b[32], const char Q[64], char R[64]) {
	u64 k1[4], k2[4];
	ec_load_k(a, k1);
	ec_load_k(b, k2);

	// Validate keys
	if (invalid_key(k1) || invalid_key(k2)) {
		return false;
	}

	// Load point
	ecpt_affine p2, r;
	ec_load_xy((const u8*)Q, p2);

	// Validate point
	if (!ec_valid(p2.x, p2.y)) {
		return false;
	}

	// Run the math routine
	ec_simul_gen(k1, k2, p2, r);

	// Save result endian-neutral
	ec_save_xy(r, (u8*)R);

	return true;
}

bool snowshoe_simul(const char a[32], const char P[64], const char b[32], const char Q[64], char R[64]) {
	u64 k1[4], k2[4];
	ec_load_k(a, k1);
	ec_load_k(b, k2);

	// Validate keys
	if (invalid_key(k1) || invalid_key(k2)) {
		return false;
	}

	// Load points
	ecpt_affine p1, p2, r;
	ec_load_xy((const u8*)P, p1);
	ec_load_xy((const u8*)Q, p2);

	// Validate points
	if (!ec_valid(p1.x, p1.y) || !ec_valid(p2.x, p2.y)) {
		return false;
	}

	// Run the math routine
	ec_simul(k1, p1, k2, p2, r);

	// Save result endian-neutral
	ec_save_xy(r, (u8*)R);

	return true;
}

#ifdef __cplusplus
}
#endif

