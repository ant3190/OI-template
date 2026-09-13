#pragma once

#include "geometry/Line.hpp"

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