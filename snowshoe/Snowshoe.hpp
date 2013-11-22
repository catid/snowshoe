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
 * R = kG
 *
 * Multiply by generator point
 */
void ec_mul_gen(const u32 k[8], ecpt &R);

/*
 * R = kP
 *
 * Multiply by variable point
 */
void ec_mul(const u32 k[8], ecpt &P, ecpt &R);

/*
 * R = aP + bQ
 *
 * Simultaneously multiply two points and return the sum
 */
void ec_simul(const u32 a[8], const ecpt &P, const u32 b[8], const ecpt &Q, ecpt &R);


#ifdef SNOWSHOE_UNIT_TESTS

void ec_mul_gen_ref1(const u32 k[8], ecpt &R);
void ec_mul_gen_ref2(const u32 k[8], ecpt &R);
void ec_mul_gen_ref3(const u32 k[8], ecpt &R);
void test_ec_mul_gen();

#endif // SNOWSHOE_UNIT_TESTS


} // namespace snowshoe

} // namespace cat

#endif // CAT_SNOWSHOE_HPP

