#pragma once

#include <algorithm>
#include <cassert>
#include <functional>
#include <utility>
#include <vector>
#include "basics/Assert.hpp"

template <class S, class Op = std::plus<>>
struct SegTree {
public:
	SegTree() : n(-1), op(), e() {}
	SegTree(int n, Op op = Op(), S e = S()) :
		n((ASSERT(n > 0), n)), ht(n == 1 ? 0 : 32 - __builtin_clz(n - 1)), m(1 << ht),
		op(std::move(op)), e(std::move(e)), tr(m << 1, this->e) {}
	SegTree(int n, S *v, Op op = Op(), S e = S()) :
		SegTree(n, std::move(op), std::move(e)) {
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
		tr[k + m - 1] = op(tr[k + m - 1], x);
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
		if (l > r) { return e; }
		S sml = e, smr = e;
		for (int lt = l + m - 1, rt = r + m; lt < rt; lt >>= 1, rt >>= 1) {
			if (lt & 1) { sml = op(sml, tr[lt++]); }
			if (rt & 1) { smr = op(tr[--rt], smr); }
		}
		return op(sml, smr);
	}
	template <class F>
	int firstright(int k, F check) {
		ASSERT(k > 0 && k <= n + 1);
		if (check(S(e))) { return k - 1; }
		if (k == n + 1) { return n + 1; }
		int kt = k + m - 1;
		S cur = e;
		do {
			while (~kt & 1) { kt >>= 1; }
			if (check(op(cur, tr[kt]))) {
				while (kt < m) {
					if (!check(op(cur, tr[kt <<= 1]))) { cur = op(cur, tr[kt++]); }
				}
				return kt - m + 1;
			}
			cur = op(cur, tr[kt]);
			++kt;
		} while (kt ^ (kt & -kt));
		return n + 1;
	}
	template <class F>
	int firstleft(int k, F check) {
		ASSERT(k >= 0 && k <= n);
		if (check(S(e))) { return k + 1; }
		if (k == 0) { return 0; }
		int kt = k + m;
		S cur = e;
		do {
			--kt;
			while (kt > 1 && kt & 1) { kt >>= 1; }
			if (check(op(tr[kt], cur))) {
				while (kt < m) {
					if (!check(op(tr[kt = kt << 1 | 1], cur))) { cur = op(tr[kt--], cur); }
				}
				return kt - m + 1;
			}
			cur = op(tr[kt], cur);
		} while (kt ^ (kt & -kt));
		return 0;
	}
	void clear() { ASSERT(n != -1); std::fill(tr.begin(), tr.end(), e); }
	
private:
	int n, ht, m;
	Op op;
	S e;
	std::vector<S> tr;
	void pushup(int p) {
		tr[p] = op(tr[p << 1], tr[p << 1 | 1]);
	}
	void up(int k) {
		for (int kt = (k + m - 1) >> 1; kt; kt >>= 1) { pushup(kt); }
	}
};
