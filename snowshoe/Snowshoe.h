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

#ifndef CAT_SNOWSHOE_HPP
#define CAT_SNOWSHOE_HPP

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Mask a provided 256-bit random number so that it is less than q
 * and can be used as a secret key.
 */
void snowshoe_secret_gen(char k[32]);

/*
 * r = (x * y + z) (mod q)
 *
 * Only one of the inputs (x) may exceed q.
 * You may pass NULL in place of z to skip the addition.
 *
 * Preconditions:
 *	0 < x
 *	0 < y,z < q
 *
 * Returns false if one of the input parameters is invalid.
 * It is important to check the return value to avoid active attacks.
 */
bool snowshoe_mul_mod_q(const char x[32], const char y[32], const char z[32], char r[32]);

/*
 * R = k*4*G
 *
 * Multiply generator point by k
 *
 * Preconditions:
 *	0 < k < q (prime order of curve)
 *
 * Returns false if one of the input parameters is invalid.
 * It is important to check the return value to avoid active attacks.
 */
bool snowshoe_mul_gen(const char k[32], const bool constant_time, char R[64]);

/*
 * R = k*4*P
 *
 * Multiply variable point by k
 *
 * Preconditions:
 * 	0 < k < q (prime order of curve)
 *
 * Returns false if one of the input parameters is invalid.
 * It is important to check the return value to avoid active attacks.
 */
bool snowshoe_mul(const char k[32], char P[64], char R[64]);

/*
 * R = a*4*P + b*4*Q
 *
 * Preconditions:
 * 	0 < a,b < q (prime order of curve)
 *
 * Simultaneously multiply two points and return the sum
 *
 * Returns false if one of the input parameters is invalid.
 * It is important to check the return value to avoid active attacks.
 */
bool snowshoe_simul(const char a[32], const char P[64], const char b[32], const char Q[64], char R[64]);

#ifdef __cplusplus
}
#endif

#endif // CAT_SNOWSHOE_HPP

