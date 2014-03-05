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
 * -(x, y) = (-x, y)
 * -(x, y, t, z) = (-x, y, -t, z)
 * t = xy/z
 *
 * In this library, the T coordinate is often split into two pieces, T0 and T1,
 * where T = T0 * T1.
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

/*
 * So about vector extensions.
 *
 * It turns out that GCC's implementation is very slow compared to Clang's.
 * On Clang there is a nice performance boost when these are used.
 * But GCC is actually twice as slow when used in the same way.
 * Therefore, I only turn on vector extensions when in Clang.
 */

#if defined(CAT_VECTOR_EXT_CLANG) && defined(CAT_WORD_64)

# define CAT_SNOWSHOE_VECTOR_OPT /* This flag is used by the rest of the code */

typedef u64 vec_ecpt_affine CAT_VECTOR_SIZE(u64, 4*2);
typedef u64 vec_ecpt CAT_VECTOR_SIZE(u64, 4*4);

#endif

static const u32 EC_D = 109;


} // namespace cat

#endif // CAT_ECPT_HPP

