#pragma once

#include "math/ModInt.hpp"

using PolyMint = ModInt<998244353>;

namespace Polynomial {
using std::vector;
using ull = unsigned long long;

constexpr int P = 998244353, Mn = 23, M = 1 << Mn;
constexpr PolyMint G = 3, invG = 332748118;

PolyMint pow_w[M], pow_iw[M], inv[M + 1];

struct initializer {
	initializer() {
		PolyMint w = G.pow(P >> Mn), iw = invG.pow(P >> Mn);
		for (int i = M >> 1; i; i >>= 1) {
			pow_w[i] = pow_iw[i] = PolyMint::raw(1);
			for (int j = 1; j < i; ++j) {
				pow_w[i + j] = pow_w[i + j - 1] * w;
				pow_iw[i + j] = pow_iw[i + j - 1] * iw;
			}
			w *= w;
			iw *= iw;
		}

		inv[0] = PolyMint::raw(1);
		for (int i = 0; i < M; ++i) {
			inv[i + 1] = inv[i] * PolyMint::raw(i + 1);
		}
		inv[M] = inv[M].inv();
		for (int i = M; i; --i) {
			PolyMint tmp = inv[i - 1];
			inv[i - 1] = inv[i] * PolyMint::raw(i);
			inv[i] *= tmp;
		}
	}
} initializer;

void NTT_DIF(vector<PolyMint> &a) {
	ASSERT(a.size() && a.size() <= M);
	ASSERT((a.size() & (a.size() - 1)) == 0);
	int len = a.size();
	for (int i = len >> 1; i; i >>= 1) {
		PolyMint *c = pow_w + i;
		for (int j = 0; j < len; j += i << 1) {
			PolyMint *p = a.data() + j, *q = p + i;
			for (int k = 0; k < i; ++k) {
				PolyMint s = p[k], t = q[k];
				p[k] = s + t;
				q[k] = (s - t) * c[k];
			}
		}
	}
}

void NTT_DIT(vector<PolyMint> &a) {
	ASSERT(a.size() && a.size() <= M);
	ASSERT((a.size() & (a.size() - 1)) == 0);
	int len = a.size();
	for (int i = 1; i < len; i <<= 1) {
		PolyMint *c = pow_iw + i;
		for (int j = 0; j < len; j += i << 1) {
			PolyMint *p = a.data() + j, *q = p + i;
			for (int k = 0; k < i; ++k) {
				PolyMint s = p[k], t = q[k] * c[k];
				p[k] = s + t;
				q[k] = s - t;
			}
		}
	}
	for (PolyMint &x : a) {
		x = x * inv[len];
	}
}

struct Poly {
	vector<PolyMint> p;

	Poly(int b = 0) { p = vector<PolyMint>(b); }
	Poly(const vector<PolyMint> &p) : p(p) {}
	Poly(const std::initializer_list<PolyMint> &p) : p(p) {}

	void clear() { p.clear(); }
	void resize(int b, PolyMint v = 0) { ASSERT(b); p.resize(b, v); }
	int size() const { return p.size(); }
	void shrink_to_fit() { p.shrink_to_fit(); }

	PolyMint &operator[](int b) { ASSERT(b < size()); return p[b]; }
	const PolyMint &operator[](int b) const { ASSERT(b < size()); return p[b]; }

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

Poly &operator+=(Poly &a, const Poly &b) {
	ASSERT(a.size() && b.size());
	if (a.size() < b.size()) {
		a.resize(b.size());
	}
	for (int i = 0; i < (int)b.size(); ++i) {
		a[i] += b[i];
	}
	return a;
}

Poly &operator-=(Poly &a, const Poly &b) {
	ASSERT(a.size() && b.size());
	if (a.size() < b.size()) { a.resize(b.size()); }
	for (int i = 0; i < (int)b.size(); ++i) { a[i] -= b[i]; }
	return a;
}

Poly &operator*=(Poly &a, Poly b) {
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

Poly &operator*=(Poly &a, PolyMint k) {
	ASSERT(a.size());
	for (int i = 0; i < (int)a.size(); ++i) {
		a[i] = a[i] * k;
	}
	return a;
}

Poly &operator<<=(Poly &a, int b) {
	ASSERT(a.size());
	a.p.insert(a.p.begin(), b, 0);
	return a;
}

Poly &operator>>=(Poly &a, int b) {
	ASSERT(a.size());
	if (b >= a.size()) {
		a.p.assign(1, 0);
	} else {
		a.p.erase(a.p.begin(), a.p.begin() + b);
	}
	return a;
}

Poly operator+(Poly a, const Poly &b) { a += b; return a; }
Poly operator-(Poly a, const Poly &b) { a -= b; return a; }
Poly operator*(Poly a, const Poly &b) { a *= b; return a; }
Poly operator*(Poly a, PolyMint k) { a *= k; return a; }
Poly operator<<(Poly a, int b) { a <<= b; return a; }
Poly operator>>(Poly a, int b) { a >>= b; return a; }

Poly Poly::ogf2egf() const {
	ASSERT(size() && size() <= M);
	Poly a(size());
	PolyMint f = PolyMint::raw(1);
	for (int i = 0; i < (int)size(); ++i) {
		a[i] = p[i] * f;
		f *= inv[i + 1];
	}
	return a;
}

Poly Poly::egf2ogf() const {
	ASSERT(size() && size() <= M);
	Poly a(size());
	PolyMint f = PolyMint::raw(1);
	for (int i = 0; i < (int)size(); ++i) {
		a[i] = p[i] * f;
		f = f * PolyMint::raw(i + 1);
	}
	return a;
}

Poly Poly::derivative() const {
	ASSERT(size());
	if (size() == 1) {
		return {0};
	}
	Poly a(size() - 1);
	for (int i = 1; i < (int)size(); ++i) {
		a[i - 1] = p[i] * PolyMint::raw(i);
	}
	return a;
}

Poly Poly::integral() const {
	ASSERT(size() && size() <= M);
	Poly a(size() + 1);
	for (int i = 0; i < (int)size(); ++i) {
		a[i + 1] = p[i] * inv[i + 1];
	}
	return a;
}

Poly Poly::square() const {
	ASSERT(size());
	int n = (size() << 1) - 1, len = 1 << (31 ^ __builtin_clz((n << 1) - 1));
	vector<PolyMint> a = p;
	a.resize(len);
	NTT_DIF(a);
	for (int i = 0; i < len; ++i) {
		a[i] *= a[i];
	}
	NTT_DIT(a);
	a.resize(n);
	return a;
}

Poly Poly::inverse() const {
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
Poly Poly::inverse(int b) const {
	ASSERT(size());
	Poly a = p;
	a.resize(b);
	return a.inverse();
}

Poly Poly::log() const {
	ASSERT(size() && p[0] == PolyMint::raw(1));
	if (size() == 1) return {0};
	Poly a = derivative() * inverse();
	a.resize(size() - 1);
	return a.integral();
}

Poly Poly::log(int b) const {
	ASSERT(size());
	Poly a = p;
	a.resize(b);
	return a.log();
}

Poly Poly::exp() const {
	ASSERT(size() && p[0] == PolyMint::raw(0));
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

Poly Poly::exp(int b) const {
	ASSERT(size());
	Poly a = p;
	a.resize(b);
	return a.exp();
}

Poly Poly::pow(long long k) const {
	ASSERT(size() && k >= 0);
	int n = size(), t = 0;
	if (!k) {
		Poly a(n);
		a[0] = 1;
		return a;
	}
	while (t < n && !p[t]) { ++t; }
	if (t == n || (t && k > (n - 1) / t)) { return Poly(n); }
	PolyMint c = p[t], ic = c.inv();
	int m = n - t * k;
	Poly a(m);
	for (int i = 0; i < m; ++i) { a[i] = p[i + t] * ic; }
	a = (a.log() * PolyMint(k)).exp() * c.pow(k);
	a <<= t * k;
	return a;
}

Poly Poly::pow(long long k, int b) const {
	ASSERT(size());
	Poly a = p;
	a.resize(b);
	return a.pow(k);
}
}

using Polynomial::Poly;