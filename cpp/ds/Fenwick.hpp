#pragma once

#include "basics/Assert.hpp"
#include "utility/Identity.hpp"

template <class T, class Op = std::plus<>, class Sub = std::minus<>, auto E = nullptr>
struct Fenwick {
public:
	Fenwick() : n(-1) {}
	Fenwick(int n) : n((ASSERT(n >= 0), n)), tr(n + 1 + (n >> 10), identity()) {}

	int size() const { return n; }

	void update(int k, T x) {
		ASSERT(n != -1 && k > 0);
		for (int kt = k; kt <= n; kt += kt & -kt) {
			int p = kt + (kt >> 10);
			tr[p] = Op{}(tr[p], x);
		}
	}
	void reset(int k) {
		ASSERT(n != -1 && k > 0);
		const T e = identity();
		for (int kt = k; kt <= n; kt += kt & -kt) { tr[kt + (kt >> 10)] = e; }
	}
	T query(int k) {
		ASSERT(k >= 0 && k <= n);
		T res = identity();
		for (int kt = k; kt; kt &= (kt - 1)) {
			res = Op{}(res, tr[kt + (kt >> 10)]);
		}
		return res;
	}
	T query(int l, int r) {
		ASSERT(n != -1);
		return (l > r ? identity() : Sub{}(query(r), query(l - 1)));
	}
	void clear() { ASSERT(n != -1); std::fill(tr.begin(), tr.end(), identity()); }
	
private:
	int n;
	std::vector<T> tr;

	static constexpr T identity() { return identity_detail::get<T, E>(); }
};
