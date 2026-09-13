#pragma once

#include "DataStructure/Fenwick.hpp"

template<class T>
struct FenwickRARS {
public:
	FenwickRARS() = default;
	FenwickRARS(int n) : tr1(n), tr2(n) {}

	void reset(int k) {
		tr1.reset(k), tr2.reset(k);
	}
	void upd(int l, int r, T x) {
		if (l <= r) { upd(l, x), upd(r + 1, -x); }
	}
	T qry(int k) {
		return tr1.qry(k) * (k + 1) - tr2.qry(k);
	}
	T qry(int l, int r) { return (l > r ? T() : qry(r) - qry(l - 1)); }
	void clear() { tr1.clear(), tr2.clear(); }

private:
	Fenwick<T> tr1, tr2;

	void upd(int k, T x) {
		tr1.upd(k, x), tr2.upd(k, x * k);
	}
};
