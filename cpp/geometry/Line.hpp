#pragma once

#include "geometry/Point.hpp"

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