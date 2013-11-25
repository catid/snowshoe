#ifndef CAT_ECPT_HPP
#define CAT_ECPT_HPP

#include "fe.hpp"

namespace cat {


/*
 * Extended Twisted Edwards Group Laws [5]
 *
 * Curve: a * u * x^2 + y^2 = d * u * x^2 * y^2 /Fp^2
 *
 * p = 2^127 - 1
 * a = -1
 * d = 109
 * u = 2 + i
 * i^2 = -1
 *
 * (0, 1) is the identity element
 * (0, -1) is of order 2
 * This codebase avoids x=0 entirely.
 *
 * -(x,y) = (-x, y)
 */

struct ecpt {
	ufe x, y, t, z;
};

struct ecpt_affine {
	ufe x, y;
};


static const u32 EC_D = 109;


} // namespace cat

#endif // CAT_ECPT_HPP

