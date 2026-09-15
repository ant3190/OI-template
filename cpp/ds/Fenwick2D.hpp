#pragma once

#include "basics/Assert.hpp"

template <class T>
struct Fenwick2D {
public:
	Fenwick2D() : n(-1), m(-1) {};
	Fenwick2D(int n, int m) : 
		n((ASSERT(n >= 0 && m >= 0), n)), m(m), tr((n + 1) * (m + 13), T()) {}

	std::pair<int, int> size() const { return {n, m}; }

	void upd(int x, int y, T val) {
		ASSERT(n != -1 && x > 0);
		ASSERT(m != -1 && y > 0);
		for (int i = x; i <= n; i += i & -i) {
			for (int j = y, s = i * (m + 13); j <= m; j += j & -j) {
				tr[s + j] += val;
			}
		}
	}
	void reset(int x, int y) {
		ASSERT(n != -1 && x > 0);
		ASSERT(m != -1 && y > 0);
		for (int i = x; i <= n; i += i & -i) {
			for (int j = y, s = i * (m + 13); j <= m; j += j & -j) { tr[s + j] = T(); }
		}
	}
	T qry(int x, int y) {
		ASSERT(x >= 0 && x <= n);
		ASSERT(y >= 0 && y <= m);
		T res = T();
		for (int i = x; i; i &= (i - 1)) {
			for (int j = y, s = i * (m + 13); j; j &= (j - 1)) { res += tr[s + j]; }
		}
		return res;
	}
	T qry(int x1, int y1, int x2, int y2) {
		ASSERT(n != -1 && m != -1);
		if (x1 > x2 || y1 > y2) return T();
		return qry(x2, y2) - qry(x1 - 1, y2) - qry(x2, y1 - 1) + qry(x1 - 1, y1 - 1);
	}
	void clear() { ASSERT(n != -1 && m != -1); std::fill(tr.begin(), tr.end(), T()); }

private:
	int n, m;
	std::vector<T> tr;
};