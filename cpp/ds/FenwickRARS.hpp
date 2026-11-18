#pragma once

#include "ds/Fenwick.hpp"

template <class T, class Op = std::plus<>, class Sub = std::minus<>, class Scale = std::multiplies<>>
struct FenwickRARS {
public:
	FenwickRARS() : sub(), scale(), e() {}
	FenwickRARS(int n, Op op = Op(), T e = T(), Sub sub = Sub(), Scale scale = Scale()) :
		sub(std::move(sub)), scale(std::move(scale)), e(std::move(e)),
		tr1(n, op, this->e), tr2(n, std::move(op), this->e) {}

	int size() const { return tr1.size(); }

	void reset(int k) {
		tr1.reset(k), tr2.reset(k);
	}
	void update(int l, int r, T x) {
		ASSERT(size() != -1);
		if (l <= r) { update(l, x), update(r + 1, sub(e, x)); }
	}
	T query(int k) {
		return sub(scale(tr1.query(k), k + 1), tr2.query(k));
	}
	T query(int l, int r) {
		ASSERT(size() != -1);
		return (l > r ? e : sub(query(r), query(l - 1)));
	}
	void clear() { ASSERT(size() != -1); tr1.clear(), tr2.clear(); }

private:
	Sub sub;
	Scale scale;
	T e;
	Fenwick<T, Op> tr1, tr2;

	void update(int k, T x) {
		tr1.update(k, x), tr2.update(k, scale(x, k));
	}
};
