#pragma once

template <class S> 
struct SegTree {
public:
	SegTree() : _n(0) {}
	explicit SegTree(int n) : 
		_n(n), ht((n == 1 ? 0 : 32 - __builtin_clz(n - 1))), _m(1 << ht), tr(_m << 1) {}
	explicit SegTree(int n, S *v) : 
		_n(n), ht(n == 1 ? 0 : 32 - __builtin_clz(n - 1)), _m(1 << ht), tr(_m << 1) {
		for (int i = 1; i <= n; ++i) { tr[i + _m - 1] = v[i]; }
		for (int i = _m - 1; i >= 1; --i) { pushup(i); }
	}
	void set(int k, S x) {
		assert(k > 0 && k <= _n);
		tr[k + _m - 1] = x;
		up(k);
	}
	void upd(int k, S x) {
		assert(k > 0 && k <= _n);
		tr[k + _m - 1] = tr[k + _m - 1] + x;
		up(k);
	}
	S qry() {
		return tr[1];
	}
	S qry(int k) {
		assert(k > 0 && k <= _n);
		return tr[k + _m - 1];
	}
	S qry(int l, int r) {
		assert(l > 0 && r <= _n);
		if (l > r) { return S(); }
		S sml = S(), smr = S();
		for (int lt = l + _m - 1, rt = r + _m; lt < rt; lt >>= 1, rt >>= 1) {
			if (lt & 1) { sml = sml + tr[lt++]; }
			if (rt & 1) { smr = tr[--rt] + smr; }
		}
		return sml + smr;
	}
	template<class F>
	int firstright(int k, F check) {
		assert(k > 0 && k <= _n + 1);
		if (check(S())) { return k - 1; }
		if (k == _n + 1) { return _n + 1; }
		int kt = k + _m - 1;
		S cur = S();
		do {
			while (~kt & 1) { kt >>= 1; }
			if (check(cur + tr[kt])) {
				while (kt < _m) {
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
		int kt = k + _m - 1;
		S cur = S();
		do {
			while (kt > 1 && kt & 1) { kt >>= 1; }
			if (check(tr[kt] + cur)) {
				while (kt < _m) {
					if (!check(tr[kt = kt << 1 | 1] + cur)) { cur = tr[kt--] + cur; }
				}
				return kt - _m + 1;
			}
			cur = tr[kt] + cur;
			--kt;
		} while (kt ^ (kt & -kt));
		return 0;
	}
	void clear() { fill(tr.begin(), tr.end(), S()); }
	
private:
	int _n, ht, _m;
	vector<S> tr;
	void pushup(int p) {
		tr[p] = tr[p << 1] + tr[p << 1 | 1];
	}
	void up(int k) {
		for (int kt = (k + _m - 1) >> 1; kt; kt >>= 1) { pushup(kt); }
	}
};