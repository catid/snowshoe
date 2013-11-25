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

#include "Platform.hpp"

/*
 * Snowshoe
 *
 * Elliptic Curve Math
 */

namespace cat {

namespace snowshoe {

/*
 * Mask a provided 256-bit number so that it is less than q
 * and can be used as a secret key
 */
void ec_mask_scalar(u64 k[4]);

/*
 * R = kG
 *
 * Multiply generator point by k
 *
 * Preconditions:
 *	0 < k < q (prime order of curve)
 */
void ec_mul_gen(const u64 k[4], ecpt_affine &R);

/*
 * R = kP
 *
 * Multiply variable point by k
 *
 * Preconditions:
 * 	0 < k < q (prime order of curve)
 * 	P is in affine (X, Y) coordinates
 */
void ec_mul(const u64 k[4], ecpt_affine &P, ecpt_affine &R);

/*
 * R = aP + bQ
 *
 * Simultaneously multiply two points and return the sum
 */
void ec_simul(const u64 a[4], const ecpt_affine &P, const u64 b[4], const ecpt_affine &Q, ecpt_affine &R);


} // namespace snowshoe

} // namespace cat

#endif // CAT_SNOWSHOE_HPP

