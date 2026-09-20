#pragma once

#include "basics/Assert.hpp"

namespace Hashing {
using ull = unsigned long long;
using u128 = __uint128_t;

static constexpr ull P = (1ULL << 61) - 1;

inline ull add(ull x, ull y) { x += y; return x >= P ? x - P : x; }
inline ull sub(ull x, ull y) { return x >= y ? x - y : x + P - y; }
inline ull mul(ull x, ull y) { u128 z = (u128)x * y; return add(z & P, z >> 61); }

inline ull qpow(ull a, ull b) {
	ull res = 1;
	while (b) {
		if (b & 1) { res = mul(res, a); }
		a = mul(a, a);
		b >>= 1;
	}
	return res;
}

inline static ull rnd(ull x) {
	x += 0x9e3779b97f4a7c15ULL;
	x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
	x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
	x ^= x >> 31;
	return x;
}

inline static ull seed = rnd(std::chrono::steady_clock::now().time_since_epoch().count());
inline static ull base = rnd(seed) % (P - 512) + 256, ibase = qpow(base, P - 2);

inline ull code(ull x) {
	x = rnd(seed + x);
	x = add((x & P), (x >> 61));
	return x ? x : 1;
}

inline std::vector<ull> pw = {1}, ipw = {1};

inline void extend(int n) {
	while ((int)pw.size() <= n) { pw.push_back(mul(pw.back(), base)); }
	while ((int)ipw.size() <= n) { ipw.push_back(mul(ipw.back(), ibase)); }
}

struct Hash {
public:
	Hash() : n(0) {}
	template <class S> Hash(const S &s) { init(s); }

	ull operator()(int l, int r) const {
		ASSERT(0 <= l && l <= r && r < n);
		extend(r - l + 1);
		return sub(h[r + 1], mul(h[l], pw[r - l + 1]));
	}
	ull operator()() const { return (*this)(0, n - 1); }

	static ull merge(ull x, ull y, int len_y) {
		extend(len_y);
		return add(mul(x, pw[len_y]), y);
	}

	template <class T>
	void append(T c) { ++n, h.push_back(add(mul(h.back(), base), code(c))); }

	void reserve(int n) { h.reserve(n + 1); }

private:
	std::vector<ull> h = {0};

	int n;

	template <class S>
	void init(const S &s) {
		n = 0;
		reserve((int)s.size());
		for (auto c : s) { append(c); }
	}
};

struct RevHash {
public:
	RevHash() : n(0) {}
	template <class S> RevHash(const S &s) { init(s); }

	ull operator()(int l, int r) const {
		ASSERT(0 <= l && l <= r && r < n);
		extend(r - l + 1);
		return mul(sub(h[r + 1], h[l]), ipw[l]);
	}
	ull operator()() const { return (*this)(0, n - 1); }

	static ull merge(ull x, ull y, int len_x) {
		extend(len_x);
		return add(x, mul(y, pw[len_x]));
	}

	template <class T>
	void append(T c) { extend(n++), h.push_back(add(h.back(), mul(code(c), pw[n - 1]))); }

	void reserve(int n) { h.reserve(n + 1); }

private:
	std::vector<ull> h = {0};

	int n;

	template <class S>
	void init(const S &s) {
		n = 0;
		reserve((int)s.size());
		for (auto c : s) { append(c); }
	}
};
}

using Hashing::Hash;
using Hashing::RevHash;