#pragma once

#include "basics/Line.hpp"

template<class P> std::vector<P> convex_hull(std::vector<P> v) {
	std::sort(v.begin(), v.end());
	v.erase(std::unique(v.begin(), v.end()), v.end());
	int n = v.size();
	if (n <= 1) { return v; }
	std::vector<P> h;
	h.reserve(n + 1);
	for (int i = 0; i < n; ++i) {
		while ((int)h.size() > 1 && side_of(h.end()[-2], h.back(), v[i]) <= 0) {
			h.pop_back();
		}
		h.push_back(v[i]);
	}
	for (int i = n - 2, t = (int)h.size(); i >= 0; --i) {
		while ((int)h.size() > t && side_of(h.end()[-2], h.back(), v[i]) <= 0) {
			h.pop_back();
		}
		h.push_back(v[i]);
	}
	h.pop_back();
	return h;
} 

template<class P> bool in_convex(const std::vector<P>& v, P p, bool incl = 1) {
	int n = v.size();
	if (n == 0) { return 0; }
	if (n == 1) { return incl && p == v[0]; }
	if (n == 2) { return incl && on_seg(v[0], v[1], p); }
	int a = side_of(v[0], v[1], p), b = side_of(v[0], v[n - 1], p);
	if (a < 0 || b > 0) { return 0; }
	P q = p - v[0];
	bool e = !a && sgn((v[1] - v[0]) * q) > 0;
	int l = 1, r = n - 1;
	while (l + 1 < r) {
		int mid = (l + r) >> 1;
		P d = v[mid] - v[0];
		int c = sgn(d ^ q);
		if (c > 0 || (!c && sgn(d * q) == 2 * e - 1)) {
			l = mid;
		} else {
			r = mid;
		}
	}
	if (!a || !b) { return incl && on_seg(v[0], v[e ? l : r], p); }
	return side_of(v[l], v[r], p) >= (int)!incl;
}

template<class P> int tangent_point(const std::vector<P> &v, P d) {
	ASSERT(!v.empty() && d != P());
	int n = v.size(), l = 0, r = n - 1;
	if (n == 1) { return 0; }
	P vl = v[1] - v[0];
	int tl = sgn(vl ^ d);
	if (tl < 0 && sgn((v[0] - v[n - 1]) ^ d) >= 0) { return 0; }
	while (l + 1 < r) {
		int mid = (l + r) >> 1;
		P vm = v[mid + 1] - v[mid];
		int tm = sgn(vm ^ d), c = sgn(vl ^ vm);
		if (c ? (tl >= 0) + (tm < 0) + (c < 0) > 1 : sgn(vl * vm) < 0 ? 
			tl > 0 || (!tl && sgn(vl * d) > 0) : sgn(vl * (v[mid] - v[l])) < 0) {
			r = mid;
		} else {
			vl = vm, tl = tm;
			l = mid;
		}
	}
	return r;
}