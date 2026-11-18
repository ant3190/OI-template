#pragma once

#include <cstddef>
#include <type_traits>

namespace identity_detail {

template <class T, auto Value, bool One = false>
constexpr T get() {
	if constexpr (std::is_same_v<decltype(Value), std::nullptr_t>) {
		if constexpr (One) { return T(1); }
		else { return T(); }
	} else if constexpr (std::is_invocable_r_v<T, decltype(Value)>) {
		return Value();
	} else {
		return T(Value);
	}
}

}
