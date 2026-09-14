#pragma once

#include "basics/Assert.hpp"
#include "math/ModInt.hpp"

template<class mint>
struct Comb {
	Comb() {}
	Comb(int n) { init(n); }

	void init(int n) {
		ASSERT(n >= 0);
		fac.resize(n + 1), ifac.resize(n + 1);
		fac[0] = 1;
		for (int i = 0; i < n; ++i) { fac[i + 1] = fac[i] * mint::raw(i + 1); }
		ifac[n] = fac[n].inv();
		for (int i = n; i > 0; --i) { ifac[i - 1] = ifac[i] * mint::raw(i); }
	}
	mint operator()(int n, int m) const {
		ASSERT(n < (int)fac.size());
		return (m < 0 || m > n ? mint(0) : fac[n] * ifac[m] * ifac[n - m]);
	}
	mint P(int n, int m) const {
		ASSERT(n < (int)fac.size());
		return (m < 0 || m > n ? mint(0) : fac[n] * ifac[n - m]);
	}

private:
	std::vector<mint> fac, ifac;
};