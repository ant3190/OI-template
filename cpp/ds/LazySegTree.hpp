#pragma once

#include <algorithm>
#include <cassert>
#include <functional>
#include <utility>
#include <vector>
#include "basics/Assert.hpp"

template <class S, class T, class Op = std::plus<>, class Apply = std::plus<>, class Compose = std::plus<>>
struct LazySegTree {
public:
	LazySegTree() : n(-1), op(), apply(), compose(), e() {}
	LazySegTree(int n, Op op = Op(), Apply apply = Apply(), Compose compose = Compose(), S e = S()) :
		n((ASSERT(n > 0), n)), ht(n == 1 ? 0 : 32 - __builtin_clz(n - 1)), m(1 << ht),
		op(std::move(op)), apply(std::move(apply)), compose(std::move(compose)), e(std::move(e)),
		tr(m << 1, this->e), tag(m), vis(m) {}
	LazySegTree(int n, S *v, Op op = Op(), Apply apply = Apply(), Compose compose = Compose(), S e = S()) :
		LazySegTree(n, std::move(op), std::move(apply), std::move(compose), std::move(e)) {
		for (int i = 1; i <= n; ++i) { tr[i + m - 1] = v[i]; }
		for (int i = m - 1; i >= 1; --i) { pushup(i); }
	}

	int size() const { return n; }

	void set(int k, S x) {
		ASSERT(k > 0 && k <= n);
		down(k);
		tr[k + m - 1] = x;
		up(k);
	}
	void update(int k, T x) {
		ASSERT(k > 0 && k <= n);
		down(k);
		pushtag(k + m - 1, x);
		up(k);
	}
	void update(int l, int r, T x) {
		ASSERT(l > 0 && r <= n);
		if (l > r) { return ; }
		down(l), down(r);
		for (int lt = l + m - 1, rt = r + m; lt < rt; lt >>= 1, rt >>= 1) {
			if (lt & 1) { pushtag(lt++, x); }
			if (rt & 1) { pushtag(--rt, x); }
		}
		up(l), up(r);
	}
	S query() {
		ASSERT(n != -1);
		return tr[1];
	}
	S query(int k) {
		ASSERT(k > 0 && k <= n);
		S res = tr[k + m - 1];
		// A pending ancestor tag is newer than any tag below it.
		for (int p = (k + m - 1) >> 1; p; p >>= 1) {
			if (vis[p]) { res = apply(res, tag[p]); }
		}
		return res;
	}
	S query(int l, int r) {
		ASSERT(l > 0 && r <= n);
		if (l > r) { return e; }
		down(l), down(r);
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
		down(k);
		int kt = k + m - 1;
		S cur = e;
		do {
			while (~kt & 1) { kt >>= 1; }
			if (check(op(cur, tr[kt]))) {
				while (kt < m) {
					pushdown(kt);
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
		down(k);
		int kt = k + m;
		S cur = e;
		do {
			--kt;
			while (kt > 1 && kt & 1) { kt >>= 1; }
			if (check(op(tr[kt], cur))) {
				while (kt < m) {
					pushdown(kt);
					if (!check(op(tr[kt = kt << 1 | 1], cur))) { cur = op(tr[kt--], cur); }
				}
				return kt - m + 1;
			}
			cur = op(tr[kt], cur);
		} while (kt ^ (kt & -kt));
		return 0;
	}
	void clear() { 
		ASSERT(n != -1);
		std::fill(tr.begin(), tr.end(), e);
		std::fill(tag.begin(), tag.end(), T());
		std::fill(vis.begin(), vis.end(), 0);
	}
	
private:
	int n, ht, m;
	Op op;
	Apply apply;
	Compose compose;
	S e;
	std::vector<S> tr;
	std::vector<T> tag;
	std::vector<int> vis;
	
	void pushtag(int p, const T &x) {
		tr[p] = apply(tr[p], x);
		if (p < m) {
			tag[p] = vis[p] ? compose(tag[p], x) : x;
			vis[p] = 1;
		}
	}
	void pushup(int p) {
		if (!vis[p]) { tr[p] = op(tr[p << 1], tr[p << 1 | 1]); }
	}
	void pushdown(int p) {
		if (vis[p]) {
			pushtag(p << 1, tag[p]), pushtag(p << 1 | 1, tag[p]);
			tag[p] = T(), vis[p] = 0;
		}
	}
	void up(int k) {
		for (int kt = (k + m - 1) >> 1; kt; kt >>= 1) { pushup(kt); }
	}
	void down(int k) {
		int kt = (k + m - 1);
		for (int i = ht; i >= 1; --i) {
			if (kt >> i) { pushdown(kt >> i); }
		}
	}
};
