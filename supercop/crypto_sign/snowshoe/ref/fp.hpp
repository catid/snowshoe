#ifndef CAT_FP_HPP
#define CAT_FP_HPP

#include "Platform.hpp"
#include "EndianNeutral.hpp"

namespace cat {


// Use builtin 128-bit type
union ufp {
	u64 i[2];
	u128 w;
};


} // namespace cat

#endif // CAT_FP_HPP

