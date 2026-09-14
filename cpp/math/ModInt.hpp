#pragma once

#include "basics/Assert.hpp"

namespace modint_internal {

constexpr long long pow_mod(long long x, long long n, int m) {
	if (m == 1) { return 0; }
	unsigned long long r = 1, y = (x % m + m) % m;
	while (n) {
		if (n & 1) { r = r * y % m; }
		y = y * y % m;
		n >>= 1;
	}
	return r;
}

constexpr bool is_prime(int n) {
	if (n <= 1) { return false; }
	if (n == 2 || n == 7 || n == 61) { return true; }
	if (!(n & 1)) { return false; }
	long long d = n - 1;
	while (!(d & 1)) { d >>= 1; }
	constexpr long long base[3] = {2, 7, 61};
	for (long long a : base) {
		long long t = d, y = pow_mod(a, t, n);
		while (t != n - 1 && y != 1 && y != n - 1) {
			y = y * y % n;
			t <<= 1;
		}
		if (y != n - 1 && !(t & 1)) { return false; }
	}
	return true;
}

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
	constexpr mint inv() const {
		if constexpr (prime) {
			ASSERT(x);
			return pow(P - 2);
		} else {
			auto [g, y] = modint_internal::inv_gcd(x, P);
			ASSERT(g == 1);
			return y;
		}
	}
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

private:
	static constexpr bool prime = modint_internal::is_prime(P);
};