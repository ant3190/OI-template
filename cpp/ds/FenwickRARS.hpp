#pragma once

#include "ds/Fenwick.hpp"

template <class T>
struct FenwickRARS {
public:
	FenwickRARS() {};
	FenwickRARS(int n) : tr1(n), tr2(n) {}

	int size() const { return tr1.size(); }

	void reset(int k) {
		tr1.reset(k), tr2.reset(k);
	}
	void update(int l, int r, T x) {
		ASSERT(size() != -1);
		if (l <= r) { update(l, x), update(r + 1, -x); }
	}
	T query(int k) {
		return tr1.query(k) * (k + 1) - tr2.query(k);
	}
	T query(int l, int r) {
		ASSERT(size() != -1);
		return (l > r ? T() : query(r) - query(l - 1));
	}
	void clear() { ASSERT(size() != -1); tr1.clear(), tr2.clear(); }

private:
	Fenwick<T> tr1, tr2;

	void update(int k, T x) {
		tr1.update(k, x), tr2.update(k, x * k);
	}
};
