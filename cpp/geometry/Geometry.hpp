#pragma once

#include "basics/Assert.hpp"

namespace Geometry {
template<class T> constexpr T eps_v = static_cast<T>(1e-9L);

template<class T>
int sgn(T x) { return (x > eps_v<T>) - (x < -eps_v<T>); }

template<class T> 
struct Point {
	using P = Point;

	T x, y;
	Point() : x(), y() {}
	Point(T x, T y) : x(x), y(y) {}

	P& operator+=(const P& p) { x += p.x, y += p.y; return *this; }
	P& operator-=(const P& p) { x -= p.x, y -= p.y; return *this; }
	bool operator<(const P& p) const { return std::tie(x, y) < std::tie(p.x, p.y); }
	bool operator==(const P& p) const { return std::tie(x, y) == std::tie(p.x, p.y); }
	bool operator!=(const P& p) const { return !(*this == p); }
	P operator+(const P& p) const { return P(x + p.x, y + p.y); }
	P operator-(const P& p) const { return P(x - p.x, y - p.y); }
	T operator^(const P& p) const { return x * p.y - y * p.x; }
	T operator*(const P& p) const { return x * p.x + y * p.y; }
	P operator*(T d) const { return P(x * d, y * d); }
	P operator/(T d) const { return P(x / d, y / d); }
	friend P operator*(T d, const P& p) { return p * d; }

	T len2() const { return x * x + y * y; }
	double len() const { return std::sqrt((double)len2()); }
	double angle() const { return std::atan2((double)y, (double)x); }
};

template<class T> T dist2(Point<T> a, Point<T> b) { return (b - a).len2(); }
template<class P> double dist(P a, P b) { return (b - a).len(); }

// 1/0/-1 : Left/On line/Right
template<class P> int side_of(P s, P e, P p) { return sgn((e - s) ^ (p - s)); }
template<class P> double line_dist(P s, P e, P p) { ASSERT(s != e); return std::fabs((double)((e - s) ^ (p - s)) / (e - s).len()); }
template<class P> double seg_dist(P s, P e, P p) {
	if (s == e) { return (p - s).len(); }
	P v1 = e - s, v2 = p - s, v3 = p - e;
	if (sgn(v1 * v2) < 0) { return v2.len(); }
	if (sgn(v1 * v3) > 0) { return v3.len(); }
	return line_dist(s, e, p);
}

template<class P> bool on_seg(P s, P e, P p, bool incl = 1) { return !side_of(s, e, p) && sgn((s - p) * (e - p)) < (int)incl; }
template<class P> bool line_seg_cross(P s1, P e1, P s2, P e2, bool incl = 1) { 
	ASSERT(s1 != e1);
	return side_of(s1, e1, s2) * side_of(s1, e1, e2) < (int)incl;
}
template<class P> bool seg_cross(P s1, P e1, P s2, P e2, bool incl = 1) {
	if (!sgn((e1 - s1) ^ (e2 - s2))) { return incl && (on_seg(s1, e1, s2) || on_seg(s1, e1, e2) || on_seg(s2, e2, s1) || on_seg(s2, e2, e1)); }
	return line_seg_cross(s1, e1, s2, e2, incl) && line_seg_cross(s2, e2, s1, e1, incl);
}

// 1/0/-1 : Intersect/Parallel/Coincident
template<class P> std::pair<int, P> line_inter(P s1, P e1, P s2, P e2) {
	ASSERT(s1 != e1 && s2 != e2);
	auto l = (e1 - s1) ^ (e2 - s2);
	if (!sgn(l)) {
		return {-(side_of(s1, e1, s2) == 0), P(0, 0)};
	}
	auto p = (e1 - s2) ^ (e2 - s2), q = (e2 - s2) ^ (s1 - s2);
	return {1, (s1 * p + e1 * q) / l};
}
template<class P> P line_proj(P s, P e, P p) {
	ASSERT(s != e);
	P v = e - s;
	return s + v * ((p - s) * v / v.len2());
}
template<class P> P line_reflect(P s, P e, P p) { return line_proj(s, e, p) * 2 - p; }

template<class T> T area2(const std::vector<Point<T>>& v) {
	int n = v.size();
	T res = T();
	for (int i = 0, j = n - 1; i < n; j = i++) {
		res += v[j] ^ v[i];
	}
	return res;
}

template<class P> bool in_polygon(const std::vector<P>& v, P p, bool incl = 1) {
	int cnt = 0, n = v.size();
	for (int i = 0, j = n - 1; i < n; j = i++) {
		if (on_seg(v[j], v[i], p)) {
			return incl;
		}
		cnt ^= ((p.y < v[i].y) - (p.y < v[j].y)) * side_of(v[j], v[i], p) > 0;
	}
	return cnt;
}

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

template<class P> void polar_sort(std::vector<P> &v, P o = P()) {
	auto half = [](const P& p) {
		return p.y < 0 || (p.y == 0 && p.x < 0);
	};
	std::sort(v.begin(), v.end(), [&](const P& a, const P& b) {
		P qa = a - o, qb = b - o;
		int ha = half(qa), hb = half(qb);
		if (ha != hb) { return ha < hb; }
		auto c = qa ^ qb;
		return c != 0 ? c > 0 : qa.len2() < qb.len2();
	});
}
}