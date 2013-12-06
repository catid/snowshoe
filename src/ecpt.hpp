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

struct ecpt_z1 {
	ufe x, y, t;
};

#ifdef CAT_HAS_VECTOR_EXTENSIONS
#define CAT_HAS_VECTOR_EXTENSIONS
#endif

#ifdef SNOWSHOE_VECTORIZE_LUT

typedef u64 vec_ecpt_affine __attribute__((ext_vector_type(4*2)));
typedef u64 vec_ecpt_z1 __attribute__((ext_vector_type(4*3)));
typedef u64 vec_ecpt __attribute__((ext_vector_type(4*4)));

#endif

static const u32 EC_D = 109;


} // namespace cat

#endif // CAT_ECPT_HPP

