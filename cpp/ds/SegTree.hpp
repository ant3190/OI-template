#pragma once

#include "basics/Assert.hpp"
#include "utility/Identity.hpp"

template <class S, class Op = std::plus<>, auto E = nullptr>
struct SegTree {
public:
	SegTree() : n(-1) {}
	SegTree(int n) :
		n((ASSERT(n > 0), n)), ht(n == 1 ? 0 : 32 - __builtin_clz(n - 1)), m(1 << ht),
		tr(m << 1, identity()) {}
	SegTree(int n, S *v) : SegTree(n) {
		for (int i = 1; i <= n; ++i) { tr[i + m - 1] = v[i]; }
		for (int i = m - 1; i >= 1; --i) { pushup(i); }
	}

	int size() const { return n; }

	void set(int k, S x) {
		ASSERT(k > 0 && k <= n);
		tr[k + m - 1] = x;
		up(k);
	}
	void update(int k, S x) {
		ASSERT(k > 0 && k <= n);
		tr[k + m - 1] = Op{}(tr[k + m - 1], x);
		up(k);
	}
	S query() {
		ASSERT(n != -1);
		return tr[1];
	}
	S query(int k) {
		ASSERT(k > 0 && k <= n);
		return tr[k + m - 1];
	}
	S query(int l, int r) {
		ASSERT(l > 0 && r <= n);
		if (l > r) { return identity(); }
		S sml = identity(), smr = identity();
		for (int lt = l + m - 1, rt = r + m; lt < rt; lt >>= 1, rt >>= 1) {
			if (lt & 1) { sml = Op{}(sml, tr[lt++]); }
			if (rt & 1) { smr = Op{}(tr[--rt], smr); }
		}
		return Op{}(sml, smr);
	}
	template <class F>
	int firstright(int k, F check) {
		ASSERT(k > 0 && k <= n + 1);
		if (check(identity())) { return k - 1; }
		if (k == n + 1) { return n + 1; }
		int kt = k + m - 1;
		S cur = identity();
		do {
			while (~kt & 1) { kt >>= 1; }
			if (check(Op{}(cur, tr[kt]))) {
				while (kt < m) {
					if (!check(Op{}(cur, tr[kt <<= 1]))) { cur = Op{}(cur, tr[kt++]); }
				}
				return kt - m + 1;
			}
			cur = Op{}(cur, tr[kt]);
			++kt;
		} while (kt ^ (kt & -kt));
		return n + 1;
	}
	template <class F>
	int firstleft(int k, F check) {
		ASSERT(k >= 0 && k <= n);
		if (check(identity())) { return k + 1; }
		if (k == 0) { return 0; }
		int kt = k + m;
		S cur = identity();
		do {
			--kt;
			while (kt > 1 && kt & 1) { kt >>= 1; }
			if (check(Op{}(tr[kt], cur))) {
				while (kt < m) {
					if (!check(Op{}(tr[kt = kt << 1 | 1], cur))) { cur = Op{}(tr[kt--], cur); }
				}
				return kt - m + 1;
			}
			cur = Op{}(tr[kt], cur);
		} while (kt ^ (kt & -kt));
		return 0;
	}
	void clear() { ASSERT(n != -1); std::fill(tr.begin(), tr.end(), identity()); }
	
private:
	int n, ht, m;
	std::vector<S> tr;

	static constexpr S identity() { return identity_detail::get<S, E>(); }

	void pushup(int p) {
		tr[p] = Op{}(tr[p << 1], tr[p << 1 | 1]);
	}
	void up(int k) {
		for (int kt = (k + m - 1) >> 1; kt; kt >>= 1) { pushup(kt); }
	}
};
