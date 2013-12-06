/*
	Copyright (c) 2013 Christopher A. Taylor.  All rights reserved.

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

#ifndef CAT_SNOWSHOE_HPP
#define CAT_SNOWSHOE_HPP

#ifdef __cplusplus
extern "C" {
#endif

#define SNOWSHOE_VERSION 6

/*
 * Verify binary compatibility with the Snowshoe API on startup.
 *
 * Example:
 * 	if (snowshoe_init()) throw "Update snowshoe static library";
 *
 * Returns 0 on success.
 * Returns non-zero if the API level does not match.
 */
int _snowshoe_init(int expected_version);
#define snowshoe_init() _snowshoe_init(SNOWSHOE_VERSION)

/*
 * Mask a provided 256-bit random number so that it is less than q
 * and can be used as a secret key.
 */
void snowshoe_secret_gen(char k[32]);

/*
 * r = (x * y + z) (mod q)
 *
 * You may pass NULL in place of z to skip the addition.
 */
void snowshoe_mul_mod_q(const char x[32], const char y[32], const char z[32], char r[32]);

/*
 * r = x (mod q)
 */
void snowshoe_mod_q(const char x[64], char r[32]);

/*
 * R = -P
 *
 * Negate the given input point and store it in R
 */
void snowshoe_neg(const char P[64], char R[64]);

/*
 * Returns 0 if the input point is valid.
 * Returns non-zero if the input point is not on the curve.
 */
int snowshoe_valid(const char P[64]);

/*
 * R =?= 4 * P
 *
 * Validates input point P.  It is only necessary to validate one of the input
 * points since the goal is to test for equality.
 *
 * This function compares the given point R with the point P, after multiplying
 * P by 4.  This is useful for signature verification to make the math work out
 * and to simultaneously validate the input point P.
 *
 * Note that to compare two points you can directly compare the bytes without
 * any special support from this library.
 *
 * Returns 0 if P is valid and R == 4 * P.
 * Returns non-zero if either P is invalid or they are not equivalent.
 */
int snowshoe_equals4(const char R[64], const char P4[64]);

/*
 * R = k*G
 *
 * Multiply generator point by k
 *
 * Validates input scalar k.
 *
 * Preconditions:
 *	0 < k < q (prime order of curve)
 *
 * Returns 0 on success.
 * Returns non-zero if one of the input parameters is invalid.
 * It is important to check the return value to avoid active attacks.
 */

int snowshoe_mul_gen(const char k[32], char R[64]);

/*
 * R = k*4*P
 *
 * Multiply variable point by k
 *
 * Validates input scalar k.  Validates input point P.
 *
 * Preconditions:
 * 	0 < k < q (prime order of curve)
 *
 * Returns 0 on success.
 * Returns non-zero if one of the input parameters is invalid.
 * It is important to check the return value to avoid active attacks.
 */
int snowshoe_mul(const char k[32], const char P[64], char R[64]);

/*
 * R = a*4*G + b*4*Q
 *
 * Validates input scalars a,b.  Validates input point Q.
 *
 * WARNING: Not constant-time.  The input parameters a,b should be public knowledge.
 * This is used mainly for signature verification where the inputs are all public.
 *
 * Preconditions:
 * 	0 < a,b < q (prime order of curve)
 *
 * Simultaneously multiply two points (one being the generator point) and return the sum
 *
 * Returns 0 on success.
 * Returns non-zero if one of the input parameters is invalid.
 * It is important to check the return value to avoid active attacks.
 */
int snowshoe_simul_gen(const char a[32], const char b[32], const char Q[64], char R[64]);

/*
 * R = a*4*P + b*4*Q
 *
 * Validates input scalars a,b.  Validates input points P,Q.
 *
 * Preconditions:
 * 	0 < a,b < q (prime order of curve)
 *
 * Simultaneously multiply two points and return the sum
 *
 * Returns 0 on success.
 * Returns non-zero if one of the input parameters is invalid.
 * It is important to check the return value to avoid active attacks.
 */
int snowshoe_simul(const char a[32], const char P[64], const char b[32], const char Q[64], char R[64]);

#ifdef __cplusplus
}
#endif

#endif // CAT_SNOWSHOE_HPP

