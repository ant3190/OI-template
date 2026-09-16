#pragma once

#include "math/ModInt.hpp"

namespace Polynomial {
using std::vector;
using ull = unsigned long long;

template <int MOD>
struct NTTInfo {
	using mint = ModInt<MOD>;

	static constexpr int P = MOD, MaxMn = 23;
	static constexpr int Mn = std::min(MaxMn, (int)std::countr_zero((unsigned)(P - 1))), M = 1 << Mn;

	static consteval mint ntt_root() {
		for (mint g = 2;; ++g) {
			if (g.pow((P - 1) >> 1) == P - 1) {
				return g.pow((P - 1) >> Mn);
			}
		}
	}

	static constexpr mint W = ntt_root(), invW = W.inv();

	static inline mint pow_w[M], pow_iw[M], inv[M + 1];

	struct initializer {
		initializer() {
			mint w = W, iw = invW;
			for (int i = M >> 1; i; i >>= 1) {
				pow_w[i] = pow_iw[i] = 1;
				for (int j = 1; j < i; ++j) {
					pow_w[i + j] = pow_w[i + j - 1] * w;
					pow_iw[i + j] = pow_iw[i + j - 1] * iw;
				}
				w *= w;
				iw *= iw;
			}
			inv[0] = 1;
			for (int i = 0; i < M; ++i) {
				inv[i + 1] = inv[i] * mint::raw(i + 1);
			}
			inv[M] = inv[M].inv();
			for (int i = M; i; --i) {
				mint tmp = inv[i - 1];
				inv[i - 1] = inv[i] * mint::raw(i);
				inv[i] *= tmp;
			}
		}
	};

	static inline initializer init;
};

template <int MOD>
void NTT_DIF(vector<ModInt<MOD>> &a) {
	using Info = NTTInfo<MOD>;
	using mint = ModInt<MOD>;
	(void)Info::init;
	ASSERT(a.size() && a.size() <= Info::M);
	ASSERT((a.size() & (a.size() - 1)) == 0);
	int len = a.size();
	for (int i = len >> 1; i; i >>= 1) {
		mint *c = Info::pow_w + i;
		for (int j = 0; j < len; j += i << 1) {
			mint *p = a.data() + j, *q = p + i;
			for (int k = 0; k < i; ++k) {
				mint s = p[k], t = q[k];
				p[k] = s + t;
				q[k] = (s - t) * c[k];
			}
		}
	}
}

template <int MOD>
void NTT_DIT(vector<ModInt<MOD>> &a) {
	using Info = NTTInfo<MOD>;
	using mint = ModInt<MOD>;
	(void)Info::init;
	ASSERT(a.size() && a.size() <= Info::M);
	ASSERT((a.size() & (a.size() - 1)) == 0);
	int len = a.size();
	for (int i = 1; i < len; i <<= 1) {
		mint *c = Info::pow_iw + i;
		for (int j = 0; j < len; j += i << 1) {
			mint *p = a.data() + j, *q = p + i;
			for (int k = 0; k < i; ++k) {
				mint s = p[k], t = q[k] * c[k];
				p[k] = s + t;
				q[k] = s - t;
			}
		}
	}
	for (mint &x : a) {
		x *= Info::inv[len];
	}
}

template <int MOD>
struct Poly {
	using mint = ModInt<MOD>;
	using Info = NTTInfo<MOD>;
	static constexpr int M = Info::M;

	vector<mint> p;

	Poly(int b = 0) { p = vector<mint>(b); }
	Poly(const vector<mint> &p) : p(p) {}
	Poly(const std::initializer_list<mint> &p) : p(p) {}

	void clear() { p.clear(); }
	void resize(int b, mint v = 0) { ASSERT(b); p.resize(b, v); }
	int size() const { return p.size(); }
	void shrink_to_fit() { p.shrink_to_fit(); }

	mint &operator[](int b) { ASSERT(b < size()); return p[b]; }
	const mint &operator[](int b) const { ASSERT(b < size()); return p[b]; }

	Poly ogf2egf() const;
	Poly egf2ogf() const;
	Poly derivative() const;
	Poly integral() const;

	Poly square() const;
	Poly inverse() const;
	Poly inverse(int b) const;
	Poly log() const;
	Poly log(int b) const;
	Poly exp() const;
	Poly exp(int b) const;
	Poly pow(long long k) const;
	Poly pow(long long k, int b) const;
};

template <int MOD>
Poly<MOD> &operator+=(Poly<MOD> &a, const Poly<MOD> &b) {
	ASSERT(a.size() && b.size());
	if (a.size() < b.size()) {
		a.resize(b.size());
	}
	for (int i = 0; i < (int)b.size(); ++i) {
		a[i] += b[i];
	}
	return a;
}

template <int MOD>
Poly<MOD> &operator-=(Poly<MOD> &a, const Poly<MOD> &b) {
	ASSERT(a.size() && b.size());
	if (a.size() < b.size()) { a.resize(b.size()); }
	for (int i = 0; i < (int)b.size(); ++i) { a[i] -= b[i]; }
	return a;
}

template <int MOD>
Poly<MOD> &operator*=(Poly<MOD> &a, Poly<MOD> b) {
	ASSERT(a.size() > 0 && b.size() > 0);
	int n = a.size() + b.size() - 1, len = 1 << (31 ^ __builtin_clz((n << 1) - 1));
	a.resize(len), b.resize(len);
	NTT_DIF(a.p), NTT_DIF(b.p);
	for (int i = 0; i < len; ++i) {
		a[i] = a[i] * b[i];
	}
	NTT_DIT(a.p);
	a.resize(n);
	return a;
}

template <int MOD>
Poly<MOD> &operator*=(Poly<MOD> &a, ModInt<MOD> k) {
	ASSERT(a.size());
	for (int i = 0; i < (int)a.size(); ++i) {
		a[i] = a[i] * k;
	}
	return a;
}

template <int MOD>
Poly<MOD> &operator<<=(Poly<MOD> &a, int b) {
	ASSERT(a.size());
	a.p.insert(a.p.begin(), b, 0);
	return a;
}

template <int MOD>
Poly<MOD> &operator>>=(Poly<MOD> &a, int b) {
	ASSERT(a.size());
	if (b >= a.size()) {
		a.p.assign(1, 0);
	} else {
		a.p.erase(a.p.begin(), a.p.begin() + b);
	}
	return a;
}

template <int MOD>
Poly<MOD> operator+(Poly<MOD> a, const Poly<MOD> &b) { a += b; return a; }

template <int MOD>
Poly<MOD> operator-(Poly<MOD> a, const Poly<MOD> &b) { a -= b; return a; }

template <int MOD>
Poly<MOD> operator*(Poly<MOD> a, const Poly<MOD> &b) { a *= b; return a; }

template <int MOD>
Poly<MOD> operator*(Poly<MOD> a, ModInt<MOD> k) { a *= k; return a; }

template <int MOD>
Poly<MOD> operator<<(Poly<MOD> a, int b) { a <<= b; return a; }

template <int MOD>
Poly<MOD> operator>>(Poly<MOD> a, int b) { a >>= b; return a; }

template <int MOD>
Poly<MOD> Poly<MOD>::ogf2egf() const {
	ASSERT(size() && size() <= M);
	Poly a(size());
	mint f = 1;
	for (int i = 0; i < (int)size(); ++i) {
		a[i] = p[i] * f;
		f *= Info::inv[i + 1];
	}
	return a;
}

template <int MOD>
Poly<MOD> Poly<MOD>::egf2ogf() const {
	ASSERT(size() && size() <= M);
	Poly a(size());
	mint f = 1;
	for (int i = 0; i < (int)size(); ++i) {
		a[i] = p[i] * f;
		f *= mint::raw(i + 1);
	}
	return a;
}

template <int MOD>
Poly<MOD> Poly<MOD>::derivative() const {
	ASSERT(size());
	if (size() == 1) {
		return {0};
	}
	Poly a(size() - 1);
	for (int i = 1; i < (int)size(); ++i) {
		a[i - 1] = p[i] * mint::raw(i);
	}
	return a;
}

template <int MOD>
Poly<MOD> Poly<MOD>::integral() const {
	ASSERT(size() && size() <= M);
	Poly a(size() + 1);
	for (int i = 0; i < (int)size(); ++i) {
		a[i + 1] = p[i] * Info::inv[i + 1];
	}
	return a;
}

template <int MOD>
Poly<MOD> Poly<MOD>::square() const {
	ASSERT(size());
	int n = (size() << 1) - 1, len = 1 << (31 ^ __builtin_clz((n << 1) - 1));
	vector<mint> a = p;
	a.resize(len);
	NTT_DIF(a);
	for (int i = 0; i < len; ++i) {
		a[i] *= a[i];
	}
	NTT_DIT(a);
	a.resize(n);
	return a;
}

template <int MOD>
Poly<MOD> Poly<MOD>::inverse() const {
	ASSERT(size() && p[0]);
	Poly a({p[0].inv()});
	for (int m = 1; m < (int)size(); m <<= 1) {
		int n = std::min(m << 1, (int)size()), len = m << 1;
		Poly b = a, c(len);
		b.resize(len);
		for (int i = 0; i < n; ++i) {
			c[i] = p[i];
		}
		NTT_DIF(b.p);
		NTT_DIF(c.p);
		for (int i = 0; i < len; ++i) {
			c[i] *= b[i];
		}
		NTT_DIT(c.p);
		for (int i = 0; i < m; ++i) {
			c[i] = 0;
		}
		NTT_DIF(c.p);
		for (int i = 0; i < len; ++i) {
			c[i] *= b[i];
		}
		NTT_DIT(c.p);
		a.resize(n);
		for (int i = m; i < n; ++i) {
			a[i] = -c[i];
		}
	}
	return a;
}

template <int MOD>
Poly<MOD> Poly<MOD>::inverse(int b) const {
	ASSERT(size());
	Poly a = p;
	a.resize(b);
	return a.inverse();
}

template <int MOD>
Poly<MOD> Poly<MOD>::log() const {
	ASSERT(size() && p[0] == 1);
	if (size() == 1) return {0};
	Poly a = derivative() * inverse();
	a.resize(size() - 1);
	return a.integral();
}

template <int MOD>
Poly<MOD> Poly<MOD>::log(int b) const {
	ASSERT(size());
	Poly a = p;
	a.resize(b);
	return a.log();
}

template <int MOD>
Poly<MOD> Poly<MOD>::exp() const {
	ASSERT(size() && p[0] == 0);
	int n = size();
	if (n == 1) {
		return {1};
	}
	Poly b({1, p[1]}), c({1}), z1, z2({1, 1});
	for (int m = 2; m < n; m <<= 1) {
		Poly y = b;
		y.resize(m << 1);
		NTT_DIF(y.p);
		z1 = z2;
		Poly z(m);
		for (int i = 0; i < m; ++i) { z[i] = y[i] * z1[i]; }
		NTT_DIT(z.p);
		for (int i = 0; i < (m >> 1); ++i) { z[i] = 0; }
		NTT_DIF(z.p);
		for (int i = 0; i < m; ++i) { z[i] = z[i] * -z1[i]; }
		NTT_DIT(z.p);
		int oldc = c.size();
		c.resize(m);
		for (int i = oldc; i < m; ++i) { c[i] = z[i]; }
		z2 = c;
		z2.resize(m << 1);
		NTT_DIF(z2.p);
		Poly x(std::min(n, m));
		for (int i = 0; i < (int)x.size(); ++i) { x[i] = p[i]; }
		x = x.derivative();
		x.resize(m);
		NTT_DIF(x.p);
		for (int i = 0; i < m; ++i) { x[i] *= y[i]; }
		NTT_DIT(x.p);
		Poly db = b.derivative();
		for (int i = 0; i < m - 1; ++i) {
			x[i] -= db[i];
		}
		x.resize(m << 1);
		for (int i = 0; i < m - 1; ++i) { x[m + i] = x[i], x[i] = 0; }
		NTT_DIF(x.p);
		for (int i = 0; i < (m << 1); ++i) { x[i] *= z2[i]; }
		NTT_DIT(x.p);
		x.resize((m << 1) - 1);
		x = x.integral();
		for (int i = m; i < std::min(n, m << 1); ++i) { x[i] += p[i]; }
		for (int i = 0; i < m; ++i) { x[i] = 0; }
		NTT_DIF(x.p);
		for (int i = 0; i < (m << 1); ++i) { x[i] *= y[i]; }
		NTT_DIT(x.p);
		int oldb = b.size();
		b.resize(std::min(n, m << 1));
		for (int i = oldb; i < (int)b.size(); ++i) { b[i] = x[i]; }
	}
	b.resize(n);
	return b;
}

template <int MOD>
Poly<MOD> Poly<MOD>::exp(int b) const {
	ASSERT(size());
	Poly a = p;
	a.resize(b);
	return a.exp();
}

template <int MOD>
Poly<MOD> Poly<MOD>::pow(long long k) const {
	ASSERT(size() && k >= 0);
	int n = size(), t = 0;
	if (!k) {
		Poly a(n);
		a[0] = 1;
		return a;
	}
	while (t < n && !p[t]) { ++t; }
	if (t == n || (t && k > (n - 1) / t)) { return Poly(n); }
	mint c = p[t], ic = c.inv();
	int m = n - t * k;
	Poly a(m);
	for (int i = 0; i < m; ++i) { a[i] = p[i + t] * ic; }
	a = (a.log() * mint(k)).exp() * c.pow(k);
	a <<= t * k;
	return a;
}

template <int MOD>
Poly<MOD> Poly<MOD>::pow(long long k, int b) const {
	ASSERT(size());
	Poly a = p;
	a.resize(b);
	return a.pow(k);
}
}

using Polynomial::Poly;