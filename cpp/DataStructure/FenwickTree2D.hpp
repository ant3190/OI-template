#pragma once

template<class T>
struct FenwickTree2D {
public:
	FenwickTree2D() : _n(0), _m(0) {}
	explicit FenwickTree2D(int n, int m) : 
		_n(n), _m(m),  
		tr((n + 1) * (m + 13), T()) {}

	void upd(int x, int y, T val) {
		assert(x > 0 && x <= _n);
		assert(y > 0 && y <= _m);
		for (int i = x; i <= _n; i += i & -i) {
			for (int j = y, s = i * (_m + 13); j <= _m; j += j & -j) {
				tr[s + j] += val;
			}
		}
	}
	void reset(int x, int y) {
		assert(x > 0 && x <= _n);
		assert(y > 0 && y <= _m);
		for (int i = x; i <= _n; i += i & -i) {
			for (int j = y, s = i * (_m + 13); j <= _m; j += j & -j) { tr[s + j] = T(); }
		}
	}
	T qry(int x, int y) {
		assert(x >= 0 && x <= _n);
		assert(y >= 0 && y <= _m);
		T res = T();
		for (int i = x; i; i &= (i - 1)) {
			for (int j = y, s = i * (_m + 13); j; j &= (j - 1)) { res += tr[s + j]; }
		}
		return res;
	}
	T qry(int x1, int y1, int x2, int y2) {
		if (x1 > x2 || y1 > y2) return T();
		return qry(x2, y2) - qry(x1 - 1, y2) - qry(x2, y1 - 1) + qry(x1 - 1, y1 - 1);
	}
	void clear() { fill(tr.begin(), tr.end(), T()); }

private:
	int _n, _m;
	vector<T> tr;
};