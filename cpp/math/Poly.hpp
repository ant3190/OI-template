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

	struct RootInfo {
		mint root[Mn + 1], iroot[Mn + 1];
		mint rate2[Mn + 1], irate2[Mn + 1];
		mint rate3[Mn + 1], irate3[Mn + 1];
		mint inv2[Mn + 1];

		consteval RootInfo() {
			root[Mn] = W, iroot[Mn] = invW;
			for (int i = Mn; i; --i) {
				root[i - 1] = root[i] * root[i];
				iroot[i - 1] = iroot[i] * iroot[i];
			}
			mint prod = 1, iprod = 1;
			for (int i = 0; i <= Mn - 2; ++i) {
				rate2[i] = root[i + 2] * prod;
				irate2[i] = iroot[i + 2] * iprod;
				prod *= iroot[i + 2];
				iprod *= root[i + 2];
			}
			prod = iprod = 1;
			for (int i = 0; i <= Mn - 3; ++i) {
				rate3[i] = root[i + 3] * prod;
				irate3[i] = iroot[i + 3] * iprod;
				prod *= iroot[i + 3];
				iprod *= root[i + 3];
			}
			inv2[0] = 1;
			mint half = mint(2).inv();
			for (int i = 1; i <= Mn; ++i) {
				inv2[i] = inv2[i - 1] * half;
			}
		}
	};

	static constexpr RootInfo info{};

	static inline vector<mint> inv{0, 1};

	static void ensure_inv(int n) {
		if ((int)inv.size() > n) { return; }
		int old = inv.size();
		inv.resize(n + 1);
		for (int i = old; i <= n; ++i) {
			inv[i] = mint::raw(P - P / i) * inv[P % i];
		}
	}
};

template <int MOD>
void NTT_DIF(vector<ModInt<MOD>> &a) {
	using Info = NTTInfo<MOD>;
	using mint = ModInt<MOD>;
	ASSERT(a.size() && a.size() <= Info::M);
	ASSERT((a.size() & (a.size() - 1)) == 0);
	constexpr auto &info = Info::info;
	int len = a.size(), h = std::countr_zero((unsigned)len), dep = 0;
	while (dep < h) {
		if (h - dep == 1) {
			int p = 1 << (h - dep - 1);
			mint rot = 1;
			for (int s = 0; s < (1 << dep); ++s) {
				int offset = s << (h - dep);
				for (int i = 0; i < p; ++i) {
					mint l = a[offset + i], r = a[offset + i + p] * rot;
					a[offset + i] = l + r;
					a[offset + i + p] = l - r;
				}
				if (s + 1 != (1 << dep)) {
					rot *= info.rate2[std::countr_zero(~(unsigned)s)];
				}
			}
			++dep;
		} else {
			int p = 1 << (h - dep - 2);
			mint rot = 1, imag = info.root[2];
			for (int s = 0; s < (1 << dep); ++s) {
				mint rot2 = rot * rot, rot3 = rot2 * rot;
				int offset = s << (h - dep);
				for (int i = 0; i < p; ++i) {
					mint a0 = a[offset + i];
					mint a1 = a[offset + i + p] * rot;
					mint a2 = a[offset + i + 2 * p] * rot2;
					mint a3 = a[offset + i + 3 * p] * rot3;
					mint a1na3imag = (a1 - a3) * imag;
					a[offset + i] = a0 + a1 + a2 + a3;
					a[offset + i + p] = a0 - a1 + a2 - a3;
					a[offset + i + 2 * p] = a0 - a2 + a1na3imag;
					a[offset + i + 3 * p] = a0 - a2 - a1na3imag;
				}
				if (s + 1 != (1 << dep)) {
					rot *= info.rate3[std::countr_zero(~(unsigned)s)];
				}
			}
			dep += 2;
		}
	}
}

template <int MOD>
void NTT_DIT(vector<ModInt<MOD>> &a) {
	using Info = NTTInfo<MOD>;
	using mint = ModInt<MOD>;
	ASSERT(a.size() && a.size() <= Info::M);
	ASSERT((a.size() & (a.size() - 1)) == 0);
	constexpr auto &info = Info::info;
	int len = a.size(), h = std::countr_zero((unsigned)len), dep = h;
	while (dep) {
		if (dep == 1) {
			int p = 1 << (h - dep);
			mint irot = 1;
			for (int s = 0; s < (1 << (dep - 1)); ++s) {
				int offset = s << (h - dep + 1);
				for (int i = 0; i < p; ++i) {
					mint l = a[offset + i], r = a[offset + i + p];
					a[offset + i] = l + r;
					a[offset + i + p] = (l - r) * irot;
				}
				if (s + 1 != (1 << (dep - 1))) {
					irot *= info.irate2[std::countr_zero(~(unsigned)s)];
				}
			}
			--dep;
		} else {
			int p = 1 << (h - dep);
			mint irot = 1, iimag = info.iroot[2];
			for (int s = 0; s < (1 << (dep - 2)); ++s) {
				mint irot2 = irot * irot, irot3 = irot2 * irot;
				int offset = s << (h - dep + 2);
				for (int i = 0; i < p; ++i) {
					mint a0 = a[offset + i];
					mint a1 = a[offset + i + p];
					mint a2 = a[offset + i + 2 * p];
					mint a3 = a[offset + i + 3 * p];
					mint a2na3iimag = (a2 - a3) * iimag;
					a[offset + i] = a0 + a1 + a2 + a3;
					a[offset + i + p] = (a0 - a1 + a2na3iimag) * irot;
					a[offset + i + 2 * p] = (a0 + a1 - a2 - a3) * irot2;
					a[offset + i + 3 * p] = (a0 - a1 - a2na3iimag) * irot3;
				}
				if (s + 1 != (1 << (dep - 2))) {
					irot *= info.irate3[std::countr_zero(~(unsigned)s)];
				}
			}
			dep -= 2;
		}
	}
	mint ilen = info.inv2[h];
	for (mint &x : a) {
		x *= ilen;
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

template <int MOD> Poly<MOD> operator+(Poly<MOD> a, const Poly<MOD> &b) { a += b; return a; }
template <int MOD> Poly<MOD> operator-(Poly<MOD> a, const Poly<MOD> &b) { a -= b; return a; }
template <int MOD> Poly<MOD> operator*(Poly<MOD> a, const Poly<MOD> &b) { a *= b; return a; }
template <int MOD> Poly<MOD> operator*(Poly<MOD> a, ModInt<MOD> k) { a *= k; return a; }
template <int MOD> Poly<MOD> operator<<(Poly<MOD> a, int b) { a <<= b; return a; }
template <int MOD> Poly<MOD> operator>>(Poly<MOD> a, int b) { a >>= b; return a; }

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
	Info::ensure_inv(size());
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
