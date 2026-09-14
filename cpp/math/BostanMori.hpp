#pragma once

#include "math/Polynomial.hpp"

PolyMint bostan_mori(Poly p, Poly q, long long n) {
	while (n) {
		Poly s = q;
		for (int i = 1; i < s.size(); i += 2) { s[i] = -s[i]; }
		Poly pt = p * s, qt = q * s;
		p.resize((pt.size() - (n & 1) + 1) / 2);
		q.resize((qt.size() + 1) / 2);
		for (int i = n & 1; i < pt.size(); i += 2) {
			p[i >> 1] = pt[i];
		}
		for (int i = 0; i < qt.size(); i += 2) {
			q[i >> 1] = qt[i];
		}
		n >>= 1;
	}
	return p[0] / q[0];
}