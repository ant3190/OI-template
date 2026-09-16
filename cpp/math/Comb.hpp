#pragma once

#include "basics/Assert.hpp"
#include "math/ModInt.hpp"

template <class mint>
struct Comb {
	Comb() {}
	Comb(int n) { init(n); }

	void init(int n) {
		ASSERT(n >= 0 && n < mint::mod());
		_fac.resize(n + 1), _ifac.resize(n + 1), _inv.resize(n + 1);
		_fac[0] = 1;
		for (int i = 0; i < n; ++i) { _fac[i + 1] = _fac[i] * mint::raw(i + 1); }
		_ifac[n] = _fac[n].inv();
		for (int i = n; i > 0; --i) { _ifac[i - 1] = _ifac[i] * mint::raw(i); }
		for (int i = 1; i <= n; ++i) { _inv[i] = _ifac[i] * _fac[i - 1]; }
	}
	mint operator()(int n, int m) const {
		if (m < 0 || m > n) { return mint(0); }
		ASSERT(n < (int)_fac.size());
		return _fac[n] * _ifac[m] * _ifac[n - m];
	}
	mint P(int n, int m) const {
		if (m < 0 || m > n) { return mint(0); }
		ASSERT(n < (int)_fac.size());
		return _fac[n] * _ifac[n - m];
	}
	mint fac(int n) const {
		ASSERT(n >= 0 && n < (int)_fac.size());
		return _fac[n];
	}
	mint ifac(int n) const {
		ASSERT(n >= 0 && n < (int)_ifac.size());
		return _ifac[n];
	}
	mint inv(int n) const {
		ASSERT(n > 0 && n < (int)_inv.size());
		return _inv[n];
	}

private:
	std::vector<mint> _fac, _ifac, _inv;
};