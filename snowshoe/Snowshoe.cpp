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

#include "ecmul.cpp"

#ifdef __cplusplus
extern "C" {
#endif

static CAT_INLINE void ec_load_k(const char k_chars[32], u64 k[4]) {
	const u64 *k_raw = reinterpret_cast<const u64 *>( k_chars );

	k[0] = getLE(k_raw[0]);
	k[1] = getLE(k_raw[1]);
	k[2] = getLE(k_raw[2]);
	k[3] = getLE(k_raw[3]);
}

void snowshoe_secret_gen(char k_chars[32]) {
	u64 *k_raw = reinterpret_cast<u64 *>( k_chars );

	u64 kq[4];
	ec_load_k(k_chars, kq);

	ec_mask_scalar(kq);

	k_raw[0] = getLE(kq[0]);
	k_raw[1] = getLE(kq[1]);
	k_raw[2] = getLE(kq[2]);
	k_raw[3] = getLE(kq[3]);
}

static bool fast_validate_k(const u64 k[4]) {
	// If k is zero return failure
	u64 zero_test = k[0] | k[1] | k[2] | k[3];
	if (zero_test) {
		return false;
	}

	// If any of high 4 bits are set, k > q
	if (k[3] > 0x0fffffffffffffffULL) {
		return false;
	}

	return true;
}

bool snowshoe_mul_gen(const char k[32], char R[64]) {
	u64 kq[4];
	ec_load_k(k, kq);

	// If k seems invalid,
	if (!fast_validate_k(kq)) {
		return false;
	}

	// Run the math routine
	ecpt_affine r;
	ec_mul_gen(kq, r);

	// Save result endian-neutral
	ec_save_xy(r, (u8*)R);

	return true;
}

bool snowshoe_mul(const char k[32], char P[64], char R[64]) {
	u64 kq[4];
	ec_load_k(k, kq);

	// If k seems invalid,
	if (!fast_validate_k(kq)) {
		return false;
	}

	// Load point
	ecpt_affine p1, r;
	ec_load_xy((const u8*)P, p1);

	// Run the math routine
	ec_mul(kq, p1, r);

	// Save result endian-neutral
	ec_save_xy(r, (u8*)R);

	return true;
}

bool snowshoe_simul(const char a[32], const char P[64], const char b[32], const char Q[64], char R[64]) {
	u64 k1[4], k2[4];
	ec_load_k(a, k1);
	ec_load_k(b, k2);

	// If k seems invalid,
	if (!fast_validate_k(k1) || !fast_validate_k(k2)) {
		return false;
	}

	// Load point
	ecpt_affine p1, p2, r;
	ec_load_xy((const u8*)P, p1);
	ec_load_xy((const u8*)Q, p1);

	// Run the math routine
	ec_simul(k1, p1, k2, p2, r);

	// Save result endian-neutral
	ec_save_xy(r, (u8*)R);

	return true;
}

#ifdef __cplusplus
}
#endif

