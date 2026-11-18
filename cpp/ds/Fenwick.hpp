#pragma once

#include <algorithm>
#include <cassert>
#include <functional>
#include <utility>
#include <vector>
#include "basics/Assert.hpp"

// op is associative and commutative, with identity e. sub(a, b) removes b from a.
template <class T, class Op = std::plus<>, class Sub = std::minus<>>
struct Fenwick {
public:
	Fenwick() : n(-1), op(), sub(), e() {}
	Fenwick(int n, Op op = Op(), T e = T(), Sub sub = Sub()) :
		n((ASSERT(n >= 0), n)), op(std::move(op)), sub(std::move(sub)), e(std::move(e)),
		tr(n + 1 + (n >> 10), this->e) {}

	int size() const { return n; }

	void update(int k, T x) {
		ASSERT(n != -1 && k > 0);
		for (int kt = k; kt <= n; kt += kt & -kt) {
			int p = kt + (kt >> 10);
			tr[p] = op(tr[p], x);
		}
	}
	void reset(int k) {
		ASSERT(n != -1 && k > 0);
		for (int kt = k; kt <= n; kt += kt & -kt) { tr[kt + (kt >> 10)] = e; }
	}
	T query(int k) {
		ASSERT(k >= 0 && k <= n);
		T res = e;
		for (int kt = k; kt; kt &= (kt - 1)) {
			res = op(res, tr[kt + (kt >> 10)]);
		}
		return res;
	}
	T query(int l, int r) {
		ASSERT(n != -1);
		return (l > r ? e : sub(query(r), query(l - 1)));
	}
	void clear() { ASSERT(n != -1); std::fill(tr.begin(), tr.end(), e); }
	
private:
	int n;
	Op op;
	Sub sub;
	T e;
	std::vector<T> tr;
};
