#pragma once

template <class S, class T> 
struct LazySegTree {
public:
	LazySegTree() : _n(0) {}
	explicit LazySegTree(int n) : 
		_n(n), ht(n == 1 ? 0 : 32 - __builtin_clz(n - 1)), _m(1 << ht), 
		tr(_m << 1), tag(_m << 1), vis(_m << 1) {}
	explicit LazySegTree(int n, S *v) : 
		_n(n), ht(n == 1 ? 0 : 32 - __builtin_clz(n - 1)), _m(1 << ht), 
		tr(_m << 1), tag(_m << 1), vis(_m << 1) {
		for (int i = 1; i <= n; ++i) { tr[i + _m - 1] = v[i], vis[i + _m - 1] = 1; }
		for (int i = n - 1; i >= 1; --i) { pushup(i); }
	}
	void set(int k, S x) {
		assert(k > 0 && k <= _n);
		down(k);
		tr[k + _m - 1] = x;
		up(k);
	}
	void upd(int k, T x) {
		assert(k > 0 && k <= _n);
		down(k);
		pushtag(k + _m - 1, x);
		up(k);
	}
	void upd(int l, int r, T x) {
		assert(l > 0 && r <= _n);
		if (l > r) { return ; }
		down(l), down(r);
		for (int lt = l + _m - 1, rt = r + _m; lt < rt; lt >>= 1, rt >>= 1) {
			if (lt & 1) { pushtag(lt++, x); }
			if (rt & 1) { pushtag(--rt, x); }
		}
		up(l), up(r);
	}
	S qry() {
		return tr[1];
	}
	S qry(int k) {
		assert(k > 0 && k <= _n);
		S res = tr[k + _m - 1];
		for (int kt = k + _m - 1; kt; kt >>= 1) {
			if (vis[kt]) { res = res + tag[kt]; }
		}
		return res;
	}
	S qry(int l, int r) {
		assert(l > 0 && r <= _n);
		if (l > r) { return S(); }
		S sml = S(), smr = S();
		int lp = (l + _m - 1) >> 1, rp = (r + _m - 1) >> 1;
		for (int lt = l + _m - 1, rt = r + _m; lt < rt; lt >>= 1, rt >>= 1, lp >>= 1, rp >>= 1) {
			if (lt & 1) { sml = sml + tr[lt++]; }
			if (rt & 1) { smr = tr[--rt] + smr; }
			if (vis[lp]) { sml = sml + tag[lp]; }
			if (vis[rp]) { smr = smr + tag[rp]; }
		}
		for (int pt = lp; pt; pt >>= 1) {
			if (vis[pt]) { sml = sml + tag[pt]; }
		}
		for (int pt = rp; pt; pt >>= 1) {
			if (vis[pt]) { smr = smr + tag[pt]; }
		}
		return sml + smr;
	}
	template<class F>
	int firstright(int k, F check) {
		assert(k > 0 && k <= _n + 1);
		if (check(S())) { return k - 1; }
		if (k == _n + 1) { return _n + 1; }
		down(k);
		int kt = k + _m - 1;
		S cur = S();
		do {
			while (~kt & 1) { kt >>= 1; }
			if (check(cur + tr[kt])) {
				while (kt < _m) {
					pushdown(kt);
					if (!check(cur + tr[kt <<= 1])) { cur = cur + tr[kt++]; }
				}
				return kt - _m + 1;
			}
			cur = cur + tr[kt];
			++kt;
		} while (kt ^ (kt & -kt));
		return _n + 1;
	}
	template<class F>
	int firstleft(int k, F check) {
		assert(k >= 0 && k <= _n);
		if (check(S())) { return k + 1; }
		if (k == 0) { return 0; }
		down(k);
		int kt = k + _m - 1;
		S cur = S();
		do {
			while (kt > 1 && kt & 1) { kt >>= 1; }
			if (check(tr[kt] + cur)) {
				while (kt < _m) {
					pushdown(kt);
					if (!check(tr[kt = kt << 1 | 1] + cur)) { cur = tr[kt--] + cur; }
				}
				return kt - _m + 1;
			}
			cur = tr[kt] + cur;
			--kt;
		} while (kt ^ (kt & -kt));
		return 0;
	}
	void clear() { 
		fill(tr.begin(), tr.end(), S());
		fill(tag.begin(), tag.end(), T());
		fill(vis.begin(), vis.end(), 0);
	}
	
private:
	int _n, ht, _m;
	vector<S> tr;
	vector<T> tag;
	vector<int> vis;
	void pushtag(int p, T x) {
		vis[p] = 1;
		tr[p] = tr[p] + x, tag[p] = tag[p] + x;
	}
	void pushup(int p) {
		if (!vis[p]) { tr[p] = tr[p << 1] + tr[p << 1 | 1]; }
	}
	void pushdown(int p) {;
		if (vis[p]) {
			pushtag(p << 1, tag[p]), pushtag(p << 1 | 1, tag[p]);
			tag[p] = T(), vis[p] = 0;
		}
	}
	void up(int k) {
		for (int kt = (k + _m - 1) >> 1; kt; kt >>= 1) { pushup(kt); }
	}
	void down(int k) {
		int kt = (k + _m - 1);
		for (int i = ht; i >= 1; --i) {
			if (kt >> i) { pushdown(kt >> i); }
		}
	}
};