#pragma once

template<class T>
struct Fenwick {
public:
	Fenwick() = default;
	Fenwick(int n) : n(n), tr(n + 1 + (n >> 10)) {}

	void upd(int k, T x) {
		assert(k > 0);
		for (int kt = k; kt <= n; kt += kt & -kt) { tr[kt + (kt >> 10)] += x; }
	}
	void reset(int k) {
		assert(k > 0);
		for (int kt = k; kt <= n; kt += kt & -kt) { tr[kt + (kt >> 10)] = T(); }
	}
	T qry(int k) {
		assert(k >= 0 && k <= n);
		T res = T();
		for (int kt = k; kt; kt &= (kt - 1)) { res += tr[kt + (kt >> 10)]; }
		return res;
	}
	T qry(int l, int r) {
		return (l > r ? T() : qry(r) - qry(l - 1));
	}
	void clear() { std::fill(tr.begin(), tr.end(), T()); }
	
private:
	int n;
	std::vector<T> tr;
};