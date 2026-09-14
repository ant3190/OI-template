#pragma once

#include "basics/Assert.hpp"

template<int P>
struct ModInt {
	using mint = ModInt;
	using ull = unsigned long long;
	int x;
	static_assert(P > 0 && P <= (1 << 30)); 
	constexpr ModInt() : x(0) {}
	template <class T>
	constexpr ModInt(T v) : x(0) {
		long long y = v % P;
		if (y < 0) { y += P; }
		x = y;
	}
	static constexpr mint raw(int x) {
		mint a;
		a.x = x;
		return a;
	}
	explicit constexpr operator int() const { return x; }
	explicit constexpr operator bool() const { return x != 0; }
	
	static constexpr int mod() { return P; }

	constexpr mint &operator+=(const mint &b) {
		x += b.x;
		if (x >= P) { x -= P; }
		return *this;
	}
	constexpr mint &operator-=(const mint &b) {
		x -= b.x;
		if (x < 0) { x += P; }
		return *this;
	}
	constexpr mint &operator*=(const mint &b) {
		x = (ull)x * b.x % P;
		return *this;
	}
	constexpr mint pow(long long n) const {
		ASSERT(n >= 0);
		mint a = *this, r = 1;
		while (n) {
			if (n & 1) r *= a;
			a *= a;
			n >>= 1;
		}
		return r;
	}
	constexpr mint inv() const { ASSERT(x); return pow(P - 2); }
	constexpr mint &operator/=(const mint &b) { return *this *= b.inv(); }
	
	constexpr mint operator-() const { return raw(x ? P - x : 0); }

	constexpr mint &operator++() { if (++x == P) { x = 0; } return *this; }
	constexpr mint &operator--() { if (x == 0) { x = P; } --x; return *this; }
	constexpr mint operator++(int) { mint res = *this; ++*this; return res; }
	constexpr mint operator--(int) { mint res = *this; --*this; return res; }

	friend constexpr mint operator+(mint a, const mint &b) { return a += b; }
	friend constexpr mint operator-(mint a, const mint &b) { return a -= b; }
	friend constexpr mint operator*(mint a, const mint &b) { return a *= b; }
	friend constexpr mint operator/(mint a, const mint &b) { return a /= b; }

	friend constexpr bool operator==(const mint &a, const mint &b) { return a.x == b.x; }
	friend constexpr bool operator!=(const mint &a, const mint &b) { return a.x != b.x; }
};