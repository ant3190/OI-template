#pragma once

template <class T, class Compare = std::less<T>>
struct Tropical {
	static_assert(std::numeric_limits<T>::is_signed);

	static T inf() {
		if constexpr (std::numeric_limits<T>::has_infinity) {
			return Compare{}(T(0), T(1)) ? std::numeric_limits<T>::infinity() : -std::numeric_limits<T>::infinity();
		} else {
			return Compare{}(T(0), T(1)) ? std::numeric_limits<T>::max() / 4 : std::numeric_limits<T>::lowest() / 4;
		}
	}

	T x;

	Tropical() : x(inf()) {}
	Tropical(T x) : x(x) {}
	static Tropical one() { return Tropical(T(0)); }

	Tropical &operator+=(const Tropical &b) {
		if (Compare{}(b.x, x)) {
			x = b.x;
		}
		return *this;
	}

	Tropical &operator*=(const Tropical &b) {
		if (x == inf() || b.x == inf()) {
			x = inf();
		} else {
			x += b.x;
		}
		return *this;
	}

	friend Tropical operator+(Tropical a, const Tropical &b) { return a += b; }
	friend Tropical operator*(Tropical a, const Tropical &b) { return a *= b; }
	friend bool operator==(const Tropical &a, const Tropical &b) { return a.x == b.x; }
	friend bool operator!=(const Tropical &a, const Tropical &b) { return !(a == b); }
};