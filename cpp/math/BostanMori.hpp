#pragma once

#include "math/Polynomial.hpp"

PolyMint bostan_mori(Poly p, Poly q, long long n) {
	ASSERT(n >= 0 && p.size() && q.size() && q[0]);
	while (n) {
		Poly s = q;
		for (int i = 1; i < s.size(); i += 2) { s[i] = -s[i]; }
		Poly pt = p * s, qt = q * s;
		int pm = (pt.size() - (n & 1) + 1) / 2, qm = (qt.size() + 1) / 2;
		if (!pm) { return 0; }
		p.resize(pm), q.resize(qm);
		for (int i = n & 1; i < pt.size(); i += 2) { p[i >> 1] = pt[i]; }
		for (int i = 0; i < qt.size(); i += 2) { q[i >> 1] = qt[i]; }
		n >>= 1;
	}
	return p[0] / q[0];
}