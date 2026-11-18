#pragma once

#include "ds/Fenwick.hpp"

template <class T, class Op = std::plus<>, class Sub = std::minus<>, class Scale = std::multiplies<>, auto E = nullptr>
struct FenwickRARS {
public:
	FenwickRARS() = default;
	FenwickRARS(int n) : tr1(n), tr2(n) {}

	int size() const { return tr1.size(); }

	void reset(int k) {
		tr1.reset(k), tr2.reset(k);
	}
	void update(int l, int r, T x) {
		ASSERT(size() != -1);
		if (l <= r) { update(l, x), update(r + 1, Sub{}(identity(), x)); }
	}
	T query(int k) {
		return Sub{}(Scale{}(tr1.query(k), k + 1), tr2.query(k));
	}
	T query(int l, int r) {
		ASSERT(size() != -1);
		return (l > r ? identity() : Sub{}(query(r), query(l - 1)));
	}
	void clear() { ASSERT(size() != -1); tr1.clear(), tr2.clear(); }

private:
	Fenwick<T, Op, Sub, E> tr1, tr2;

	static constexpr T identity() { return identity_detail::get<T, E>(); }

	void update(int k, T x) {
		tr1.update(k, x), tr2.update(k, Scale{}(x, k));
	}
};
