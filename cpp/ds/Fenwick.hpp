#pragma once

#include "basics/Assert.hpp"

template <class T>
struct Fenwick {
public:
	Fenwick() : n(-1) {};
	Fenwick(int n) : n((ASSERT(n >= 0), n)), tr(n + 1 + (n >> 10)) {}

	int size() const { return n; }

	void upd(int k, T x) {
		ASSERT(n != -1 && k > 0);
		for (int kt = k; kt <= n; kt += kt & -kt) { tr[kt + (kt >> 10)] += x; }
	}
	void reset(int k) {
		ASSERT(n != -1 && k > 0);
		for (int kt = k; kt <= n; kt += kt & -kt) { tr[kt + (kt >> 10)] = T(); }
	}
	T qry(int k) {
		ASSERT(k >= 0 && k <= n);
		T res = T();
		for (int kt = k; kt; kt &= (kt - 1)) { res += tr[kt + (kt >> 10)]; }
		return res;
	}
	T qry(int l, int r) {
		ASSERT(n != -1);
		return (l > r ? T() : qry(r) - qry(l - 1));
	}
	void clear() { ASSERT(n != -1); std::fill(tr.begin(), tr.end(), T()); }
	
private:
	int n;
	std::vector<T> tr;
};