#pragma once

#include "basics/Assert.hpp"
#include "utility/Identity.hpp"

template <class S, class T, class Op = std::plus<>, class Apply = std::plus<>, class Compose = std::plus<>,
		auto E = nullptr, auto TagId = nullptr>
struct LazySegTree {
public:
	LazySegTree() : n(-1) {}
	LazySegTree(int n) { init(n); }
	LazySegTree(int n, S *v) { init(n, v); }

	void init(int N) {
		ASSERT(N > 0);
		n = N, ht = (n == 1 ? 0 : 32 - __builtin_clz(n - 1)), m = (1 << ht);
		tr = std::vector<T>(m << 1, identity()), tag = std::vector<T>(m, tag_identity()), vis = std::vector<int>(m);
	}
	void init(int n, S *v) {
		init(n);
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
			if (vis[p]) { res = Apply{}(res, tag[p]); }
		}
		return res;
	}
	S query(int l, int r) {
		ASSERT(l > 0 && r <= n);
		if (l > r) { return identity(); }
		down(l), down(r);
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
		down(k);
		int kt = k + m - 1;
		S cur = identity();
		do {
			while (~kt & 1) { kt >>= 1; }
			if (check(Op{}(cur, tr[kt]))) {
				while (kt < m) {
					pushdown(kt);
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
		down(k);
		int kt = k + m;
		S cur = identity();
		do {
			--kt;
			while (kt > 1 && kt & 1) { kt >>= 1; }
			if (check(Op{}(tr[kt], cur))) {
				while (kt < m) {
					pushdown(kt);
					if (!check(Op{}(tr[kt = kt << 1 | 1], cur))) { cur = Op{}(tr[kt--], cur); }
				}
				return kt - m + 1;
			}
			cur = Op{}(tr[kt], cur);
		} while (kt ^ (kt & -kt));
		return 0;
	}
	void clear() { 
		ASSERT(n != -1);
		std::fill(tr.begin(), tr.end(), identity());
		std::fill(tag.begin(), tag.end(), tag_identity());
		std::fill(vis.begin(), vis.end(), 0);
	}
	
private:
	int n, ht, m;
	std::vector<S> tr;
	std::vector<T> tag;
	std::vector<int> vis;
	
	static constexpr S identity() { return identity_detail::get<S, E>(); }
	static constexpr T tag_identity() { return identity_detail::get<T, TagId>(); }

	void pushtag(int p, const T &x) {
		tr[p] = Apply{}(tr[p], x);
		if (p < m) {
			tag[p] = vis[p] ? Compose{}(tag[p], x) : x;
			vis[p] = 1;
		}
	}
	void pushup(int p) {
		if (!vis[p]) { tr[p] = Op{}(tr[p << 1], tr[p << 1 | 1]); }
	}
	void pushdown(int p) {
		if (vis[p]) {
			pushtag(p << 1, tag[p]), pushtag(p << 1 | 1, tag[p]);
			tag[p] = tag_identity(), vis[p] = 0;
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
