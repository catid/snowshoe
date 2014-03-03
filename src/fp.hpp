#ifndef CAT_FP_HPP
#define CAT_FP_HPP

#include "Platform.hpp"
#include "EndianNeutral.hpp"
#include "BigMath.hpp" // 128-bit math

namespace cat {


union ufp {
	u64 i[2];
	u128 w;
};


} // namespace cat

#endif // CAT_FP_HPP

