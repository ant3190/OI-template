#pragma once

#include "basics/Assert.hpp"

template <class S> 
struct SegTree {
public:
	SegTree() = default;
	SegTree(int n) : 
		n((ASSERT(n > 0), n)), ht((n == 1 ? 0 : 32 - __builtin_clz(n - 1))), m(1 << ht), tr(m << 1) {}
	SegTree(int n, S *v) : 
		n((ASSERT(n > 0), n)), ht(n == 1 ? 0 : 32 - __builtin_clz(n - 1)), m(1 << ht), tr(m << 1) {
		for (int i = 1; i <= n; ++i) { tr[i + m - 1] = v[i]; }
		for (int i = m - 1; i >= 1; --i) { pushup(i); }
	}

	int size() const { return n; }

	void set(int k, S x) {
		ASSERT(k > 0 && k <= n);
		tr[k + m - 1] = x;
		up(k);
	}
	void upd(int k, S x) {
		ASSERT(k > 0 && k <= n);
		tr[k + m - 1] = tr[k + m - 1] + x;
		up(k);
	}
	S qry() {
		ASSERT(n != -1);
		return tr[1];
	}
	S qry(int k) {
		ASSERT(k > 0 && k <= n);
		return tr[k + m - 1];
	}
	S qry(int l, int r) {
		ASSERT(l > 0 && r <= n);
		if (l > r) { return S(); }
		S sml = S(), smr = S();
		for (int lt = l + m - 1, rt = r + m; lt < rt; lt >>= 1, rt >>= 1) {
			if (lt & 1) { sml = sml + tr[lt++]; }
			if (rt & 1) { smr = tr[--rt] + smr; }
		}
		return sml + smr;
	}
	template<class F>
	int firstright(int k, F check) {
		ASSERT(k > 0 && k <= n + 1);
		if (check(S())) { return k - 1; }
		if (k == n + 1) { return n + 1; }
		int kt = k + m - 1;
		S cur = S();
		do {
			while (~kt & 1) { kt >>= 1; }
			if (check(cur + tr[kt])) {
				while (kt < m) {
					if (!check(cur + tr[kt <<= 1])) { cur = cur + tr[kt++]; }
				}
				return kt - m + 1;
			}
			cur = cur + tr[kt];
			++kt;
		} while (kt ^ (kt & -kt));
		return n + 1;
	}
	template<class F>
	int firstleft(int k, F check) {
		ASSERT(k >= 0 && k <= n);
		if (check(S())) { return k + 1; }
		if (k == 0) { return 0; }
		int kt = k + m - 1;
		S cur = S();
		do {
			while (kt > 1 && kt & 1) { kt >>= 1; }
			if (check(tr[kt] + cur)) {
				while (kt < m) {
					if (!check(tr[kt = kt << 1 | 1] + cur)) { cur = tr[kt--] + cur; }
				}
				return kt - m + 1;
			}
			cur = tr[kt] + cur;
			--kt;
		} while (kt ^ (kt & -kt));
		return 0;
	}
	void clear() { ASSERT(n != -1); std::fill(tr.begin(), tr.end(), S()); }
	
private:
	int n = -1, ht, m;
	std::vector<S> tr;
	void pushup(int p) {
		tr[p] = tr[p << 1] + tr[p << 1 | 1];
	}
	void up(int k) {
		for (int kt = (k + m - 1) >> 1; kt; kt >>= 1) { pushup(kt); }
	}
};