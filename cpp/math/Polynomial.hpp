#pragma once

#include "basics/Polynomial.hpp"

namespace Polynomial {
using std::vector;
using ull = unsigned long long;

constexpr int mod = 998244353, G = 3, invG = 332748118, Mn = 23, M = 1 << Mn;

int pow_w[M], pow_iw[M], inv[M + 1];

template <typename _Tp>
void addmod(_Tp &x) {
	if (x >= mod) {
		x -= mod;
	}
}

int qpow(int a, int b) {
	int res = 1;
	while (b) {
		if (b & 1) {
			res = (ull)res * a % mod;
		}
		a = (ull)a * a % mod;
		b >>= 1;
	}
	return res;
}

int qinv(int a) {
	return qpow(a, mod - 2);
}

struct initializer {
	initializer() {
		int w = qpow(G, mod >> Mn), iw = qpow(invG, mod >> Mn);
		for (int i = M >> 1; i; i >>= 1) {
			pow_w[i] = pow_iw[i] = 1;
			for (int j = 1; j < i; ++j) {
				pow_w[i + j] = (ull)pow_w[i + j - 1] * w % mod;
				pow_iw[i + j] = (ull)pow_iw[i + j - 1] * iw % mod;
			}
			w = (ull)w * w % mod;
			iw = (ull)iw * iw % mod;
		}

		inv[0] = 1;
		for (int i = 0; i < M; ++i) {
			inv[i + 1] = (ull)inv[i] * (i + 1) % mod;
		}
		inv[M] = qinv(inv[M]);
		for (int i = M; i; --i) {
			int tmp = inv[i - 1];
			inv[i - 1] = (ull)inv[i] * i % mod;
			inv[i] = (ull)inv[i] * tmp % mod;
		}
	}
} initializer;

void NTT_DIF(vector<int> &a) {
	ASSERT(a.size() && a.size() <= M);
	ASSERT((a.size() & (a.size() - 1)) == 0);
	int len = a.size();
	for (int i = len >> 1; i; i >>= 1) {
		int *c = pow_w + i;
		for (int j = 0; j < len; j += i << 1) {
			int *p = a.data() + j, *q = p + i;
			for (int k = 0; k < i; ++k) {
				int s = p[k], t = q[k];
				addmod(p[k] = s + t);
				q[k] = (ull)(s + mod - t) * c[k] % mod;
			}
		}
	}
}

void NTT_DIT(vector<int> &a) {
	ASSERT(a.size() && a.size() <= M);
	ASSERT((a.size() & (a.size() - 1)) == 0);
	int len = a.size();
	for (int i = 1; i < len; i <<= 1) {
		int *c = pow_iw + i;
		for (int j = 0; j < len; j += i << 1) {
			int *p = a.data() + j, *q = p + i;
			for (int k = 0; k < i; ++k) {
				int s = p[k], t = (ull)q[k] * c[k] % mod;
				addmod(p[k] = s + t);
				addmod(q[k] = s + mod - t);
			}
		}
	}
	for (int &x : a) {
		x = (ull)x * inv[len] % mod;
	}
}

struct Poly {
	vector<int> p;

	Poly(size_t b = 0) { p = vector<int>(b); }
	Poly(const vector<int> &p) : p(p) {}
	Poly(const std::initializer_list<int> &p) : p(p) {}

	void clear() { p.clear(); }
	void resize(size_t b, int v = 0) { ASSERT(b && 0 <= v && v < mod); p.resize(b, v); }
	size_t size() const { return p.size(); }
	void shrink_to_fit() { p.shrink_to_fit(); }

	int &operator[](size_t b) { ASSERT(b < size()); return p[b]; }
	const int &operator[](size_t b) const { ASSERT(b < size()); return p[b]; }

	void print(const std::string &sep = " ", const std::string &end = "\n") const {
		for (int i = 0; i < (int)size(); ++i) {
			if (i) {
				std::cout << sep;
			}
			std::cout << p[i];
		}
		std::cout << end;
	}

	Poly ogf2egf() const;
	Poly egf2ogf() const;
	Poly derivative() const;
	Poly integral() const;

	Poly square() const;
	Poly inverse() const;
	Poly inverse(size_t b) const;
	Poly log() const;
	Poly log(size_t b) const;
	Poly exp() const;
	Poly exp(size_t b) const;
	Poly pow(int k) const;
	Poly pow(int k, size_t b) const;
};

Poly &operator+=(Poly &a, const Poly &b) {
	ASSERT(a.size() && b.size());
	if (a.size() < b.size()) {
		a.resize(b.size());
	}
	for (int i = 0; i < (int)b.size(); ++i) {
		addmod(a[i] += b[i]);
	}
	return a;
}

Poly &operator-=(Poly &a, const Poly &b) {
	ASSERT(a.size() && b.size());
	if (a.size() < b.size()) {
		a.resize(b.size());
	}
	for (int i = 0; i < (int)b.size(); ++i) {
		addmod(a[i] += mod - b[i]);
	}
	return a;
}

Poly &operator*=(Poly &a, Poly b) {
	ASSERT(a.size() > 0 && b.size() > 0);
	int n = a.size() + b.size() - 1, len = 1 << (31 ^ __builtin_clz((n << 1) - 1));
	a.resize(len), b.resize(len);
	NTT_DIF(a.p), NTT_DIF(b.p);
	for (int i = 0; i < len; ++i) {
		a[i] = (ull)a[i] * b[i] % mod;
	}
	NTT_DIT(a.p);
	a.resize(n);
	return a;
}

Poly &operator*=(Poly &a, int k) {
	ASSERT(a.size() && 0 <= k && k < mod);
	for (int i = 0; i < (int)a.size(); ++i) {
		a[i] = (ull)a[i] * k % mod;
	}
	return a;
}

Poly &operator<<=(Poly &a, size_t b) {
	ASSERT(a.size());
	a.p.insert(a.p.begin(), b, 0);
	return a;
}

Poly operator+(Poly a, const Poly &b) { a += b; return a; }
Poly operator-(Poly a, const Poly &b) { a -= b; return a; }
Poly operator*(Poly a, const Poly &b) { a *= b; return a; }
Poly operator*(Poly a, int k) { a *= k; return a; }
Poly operator<<(Poly a, size_t b) { a <<= b; return a; }

Poly &operator>>=(Poly &a, size_t b) {
	ASSERT(a.size());
	if (b >= a.size()) {
		a.p.assign(1, 0);
	} else {
		a.p.erase(a.p.begin(), a.p.begin() + b);
	}
	return a;
}

Poly operator>>(Poly a, size_t b) { a >>= b; return a; }

Poly Poly::ogf2egf() const {
	ASSERT(size() && size() <= M);
	Poly a(size());
	for (int i = 0, f = 1; i < (int)size(); ++i) {
		a[i] = (ull)p[i] * f % mod;
		f = (ull)f * inv[i + 1] % mod;
	}
	return a;
}

Poly Poly::egf2ogf() const {
	ASSERT(size() && size() <= M);
	Poly a(size());
	for (int i = 0, f = 1; i < (int)size(); ++i) {
		a[i] = (ull)p[i] * f % mod;
		f = (ull)f * (i + 1) % mod;
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
		a[i - 1] = (ull)p[i] * i % mod;
	}
	return a;
}

Poly Poly::integral() const {
	ASSERT(size() && size() <= M);
	Poly a(size() + 1);
	for (int i = 0; i < (int)size(); ++i) {
		a[i + 1] = (ull)p[i] * inv[i + 1] % mod;
	}
	return a;
}

Poly Poly::square() const {
	ASSERT(size());
	int n = (size() << 1) - 1, len = 1 << (31 ^ __builtin_clz((n << 1) - 1));
	vector<int> a = p;
	a.resize(len);
	NTT_DIF(a);
	for (int i = 0; i < len; ++i) {
		a[i] = (ull)a[i] * a[i] % mod;
	}
	NTT_DIT(a);
	a.resize(n);
	return a;
}

Poly Poly::inverse() const {
	ASSERT(size() && p[0]);
	Poly a({qinv(p[0])});
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
			c[i] = (ull)c[i] * b[i] % mod;
		}
		NTT_DIT(c.p);
		for (int i = 0; i < m; ++i) {
			c[i] = 0;
		}
		NTT_DIF(c.p);
		for (int i = 0; i < len; ++i) {
			c[i] = (ull)c[i] * b[i] % mod;
		}
		NTT_DIT(c.p);
		a.resize(n);
		for (int i = m; i < n; ++i) {
			a[i] = c[i] ? mod - c[i] : 0;
		}
	}
	return a;
}
Poly Poly::inverse(size_t b) const {
	ASSERT(size());
	Poly a = p;
	a.resize(b);
	return a.inverse();
}

Poly Poly::log() const {
	ASSERT(size() && p[0] == 1);
	if (size() == 1) return {0};
	Poly a = derivative() * inverse();
	a.resize(size() - 1);
	return a.integral();
}

Poly Poly::log(size_t b) const {
	ASSERT(size());
	Poly a = p;
	a.resize(b);
	return a.log();
}

Poly Poly::exp() const {
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
		for (int i = 0; i < m; ++i) { z[i] = (ull)y[i] * z1[i] % mod; }
		NTT_DIT(z.p);
		for (int i = 0; i < (m >> 1); ++i) { z[i] = 0; }
		NTT_DIF(z.p);
		for (int i = 0; i < m; ++i) { z[i] = (ull)z[i] * (z1[i] ? mod - z1[i] : 0) % mod; }
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
		for (int i = 0; i < m; ++i) { x[i] = (ull)x[i] * y[i] % mod; }
		NTT_DIT(x.p);
		Poly db = b.derivative();
		for (int i = 0; i < m - 1; ++i) {
			x[i] -= db[i];
			if (x[i] < 0) { x[i] += mod; }
		}
		x.resize(m << 1);
		for (int i = 0; i < m - 1; ++i) { x[m + i] = x[i], x[i] = 0; }
		NTT_DIF(x.p);
		for (int i = 0; i < (m << 1); ++i) { x[i] = (ull)x[i] * z2[i] % mod; }
		NTT_DIT(x.p);
		x.resize((m << 1) - 1);
		x = x.integral();
		for (int i = m; i < std::min(n, m << 1); ++i) { addmod(x[i] += p[i]); }
		for (int i = 0; i < m; ++i) { x[i] = 0; }
		NTT_DIF(x.p);
		for (int i = 0; i < (m << 1); ++i) { x[i] = (ull)x[i] * y[i] % mod; }
		NTT_DIT(x.p);
		int oldb = b.size();
		b.resize(std::min(n, m << 1));
		for (int i = oldb; i < (int)b.size(); ++i) { b[i] = x[i]; }
	}
	b.resize(n);
	return b;
}

Poly Poly::exp(size_t b) const {
	ASSERT(size());
	Poly a = p;
	a.resize(b);
	return a.exp();
}

Poly Poly::pow(int k) const {
	ASSERT(size() && k >= 0);
	int n = size(), t = 0;
	if (!k) {
		Poly a(n);
		a[0] = 1;
		return a;
	}
	while (t < n && !p[t]) { ++t; }
	if (t == n || (ull)t * k >= n) { return Poly(n); }
	int c = p[t], ic = qinv(c), m = n - (ull)t * k;
	Poly a(m);
	for (int i = 0; i < m; ++i) { a[i] = (ull)p[i + t] * ic % mod; }
	a = (a.log() * (k % mod)).exp() * qpow(c, k);
	a <<= (ull)t * k;
	return a;
}

Poly Poly::pow(int k, size_t b) const {
	ASSERT(size());
	Poly a = p;
	a.resize(b);
	return a.pow(k);
}
}

using Polynomial::Poly;