#pragma once

#include "basics/Assert.hpp"

namespace dynamic_modint_internal {
struct Barrett {
	using ull = unsigned long long;
	unsigned int m;
	ull im;

	Barrett(unsigned int m = 998244353) : m(m), im(ull(-1) / m + 1) {}

	unsigned int mul(unsigned int a, unsigned int b) const {
		ull z = (ull)a * b;
		ull x = (ull)((__uint128_t)z * im >> 64);
		ull y = x * m;
		return z < y ? z + m - y : z - y;
	}
};

constexpr std::pair<long long, long long> inv_gcd(long long a, long long b) {
	a %= b;
	if (a < 0) { a += b; }
	if (!a) { return {b, 0}; }
	long long s = b, t = a, m0 = 0, m1 = 1;
	while (t) {
		long long u = s / t;
		s -= t * u, m0 -= m1 * u;
		std::swap(s, t), std::swap(m0, m1);
	}
	if (m0 < 0) { m0 += b / s; }
	return {s, m0};
}
}

template <int id = 0>
struct DynModInt {
	using mint = DynModInt;
	int x;

	DynModInt() : x(0) {}
	template <class T>
	DynModInt(T v) : x(0) {
		long long y = v % mod();
		if (y < 0) { y += mod(); }
		x = y;
	}
	static mint raw(int x) {
		mint a;
		a.x = x;
		return a;
	}
	explicit operator int() const { return x; }
	explicit operator bool() const { return x != 0; }

	static int mod() { return bt.m; }
	static void set_mod(int p) {
		ASSERT(p > 0 && p <= (1 << 30));
		bt = dynamic_modint_internal::Barrett(p);
	}

	mint &operator+=(const mint &b) {
		x += b.x;
		if (x >= mod()) { x -= mod(); }
		return *this;
	}
	mint &operator-=(const mint &b) {
		x -= b.x;
		if (x < 0) { x += mod(); }
		return *this;
	}
	mint &operator*=(const mint &b) {
		x = bt.mul(x, b.x);
		return *this;
	}
	mint pow(long long n) const {
		ASSERT(n >= 0);
		mint a = *this, r = 1;
		while (n) {
			if (n & 1) { r *= a; }
			a *= a;
			n >>= 1;
		}
		return r;
	}
	mint inv() const {
		auto [g, y] = dynamic_modint_internal::inv_gcd(x, mod());
		ASSERT(g == 1);
		return y;
	}
	mint &operator/=(const mint &b) { return *this *= b.inv(); }

	mint operator-() const { return raw(x ? mod() - x : 0); }

	mint &operator++() { if (++x == mod()) { x = 0; } return *this; }
	mint &operator--() { if (x == 0) { x = mod(); } --x; return *this; }
	mint operator++(int) { mint res = *this; ++*this; return res; }
	mint operator--(int) { mint res = *this; --*this; return res; }

	friend mint operator+(mint a, const mint &b) { return a += b; }
	friend mint operator-(mint a, const mint &b) { return a -= b; }
	friend mint operator*(mint a, const mint &b) { return a *= b; }
	friend mint operator/(mint a, const mint &b) { return a /= b; }

	friend bool operator==(const mint &a, const mint &b) { return a.x == b.x; }
	friend bool operator!=(const mint &a, const mint &b) { return a.x != b.x; }

private:
	inline static dynamic_modint_internal::Barrett bt;
};
