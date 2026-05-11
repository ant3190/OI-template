#include <bits/stdc++.h>
#include <cassert>
#define INF 0x3f3f3f3f
#define INF_LL 0x3f3f3f3f3f3f3f3f
using namespace std;

using ll = long long;
using ull = unsigned long long;
using i128 = __int128_t;

mt19937 rd(chrono::steady_clock().now().time_since_epoch().count());
mt19937_64 rdll(chrono::steady_clock().now().time_since_epoch().count());

/*========== Begin ==========*/

template<class T>
struct FenwickTree {
public:
	FenwickTree() : _n(0) {}
	explicit FenwickTree(int n) : _n(n), tr(n + 1 + (n >> 10)) {}

	void upd(int k, T x) {
		assert(k > 0 && k <= _n);
		for (int kt = k; kt <= _n; kt += kt & -kt) { tr[kt + (kt >> 10)] += x; }
	}
	void reset(int k) {
		assert(k > 0 && k <= _n);
		for (int kt = k; kt <= _n; kt += kt & -kt) { tr[kt + (kt >> 10)] = T(); }
	}
	T qry(int k) {
		assert(k >= 0 && k <= _n);
		T res = T();
		for (int kt = k; kt; kt &= (kt - 1)) { res += tr[kt + (kt >> 10)]; }
		return res;
	}
	T qry(int l, int r) {
		return (l > r ? T() : qry(r) - qry(l - 1));
	}
	void clear() { fill(tr.begin(), tr.end(), T()); }
	
private:
	int _n;
	vector<T> tr;
};

template<class T>
struct FenwickTreeRARS {
public:
	FenwickTreeRARS() : _n(0) {}
	explicit FenwickTreeRARS(int n) : _n(n), tr1(n + 1 + (n >> 10)), tr2(n + 1 + (n >> 10)) {}

	void reset(int k) {
		assert(k > 0 && k <= _n);
		for (int kt = k; kt <= _n; kt += kt & -kt) { tr1[kt + (kt >> 10)] = T(); }
		for (int kt = k; kt <= _n; kt += kt & -kt) { tr2[kt + (kt >> 10)] = T(); }
	}
	void upd(int l, int r, T x) {
		if (l <= r) { upd(l, x), upd(r + 1, -x); }
	}
	T qry(int k) {
		assert(k >= 0 && k <= _n);
		T r1 = T(), r2 = T();
		for (int kt = k; kt; kt &= (kt - 1)) { r1 += tr1[kt + (kt >> 10)]; }
		for (int kt = k; kt; kt &= (kt - 1)) { r2 += tr2[kt + (kt >> 10)]; }
		return r1 * (k + 1) - r2;
	}
	T qry(int l, int r) {
		return (l > r ? T() : qry(r) - qry(l - 1));
	}
	void clear() { fill(tr1.begin(), tr1.end(), T()), fill(tr2.begin(), tr2.end(), T()); }

private:
	int _n;
	vector<T> tr1, tr2;

	void upd(int k, T x) {
		assert(k > 0 && k <= _n + 1);
		for (int kt = k; kt <= _n; kt += kt & -kt) { tr1[kt + (kt >> 10)] += x; }
		x = x * k;
		for (int kt = k; kt <= _n; kt += kt & -kt) { tr2[kt + (kt >> 10)] += x; }
	}
};

template<class T>
struct FenwickTree2D {
public:
	FenwickTree2D() : _n(0), _m(0) {}
	explicit FenwickTree2D(int n, int m) : 
		_n(n), _m(m),  
		tr((n + 1) * (m + 13), T()) {}

	void upd(int x, int y, T val) {
		assert(x > 0 && x <= _n);
		assert(y > 0 && y <= _m);
		for (int i = x; i <= _n; i += i & -i) {
			for (int j = y, s = i * (_m + 13); j <= _m; j += j & -j) {
				tr[s + j] += val;
			}
		}
	}
	void reset(int x, int y) {
		assert(x > 0 && x <= _n);
		assert(y > 0 && y <= _m);
		for (int i = x; i <= _n; i += i & -i) {
			for (int j = y, s = i * (_m + 13); j <= _m; j += j & -j) { tr[s + j] = T(); }
		}
	}
	T qry(int x, int y) {
		assert(x >= 0 && x <= _n);
		assert(y >= 0 && y <= _m);
		T res = T();
		for (int i = x; i; i &= (i - 1)) {
			for (int j = y, s = i * (_m + 13); j; j &= (j - 1)) { res += tr[s + j]; }
		}
		return res;
	}
	T qry(int x1, int y1, int x2, int y2) {
		if (x1 > x2 || y1 > y2) return T();
		return qry(x2, y2) - qry(x1 - 1, y2) - qry(x2, y1 - 1) + qry(x1 - 1, y1 - 1);
	}
	void clear() { fill(tr.begin(), tr.end(), T()); }

private:
	int _n, _m;
	vector<T> tr;
};

/*========== End ==========*/

int n, m;

int main() {
#ifdef LOCAL
	assert(freopen("test.in", "r", stdin));
	assert(freopen("test.out", "w", stdout));
#endif

	cin.tie(0)->sync_with_stdio(0);

	cin >> n >> m;
	FenwickTreeRARS<ll> sum(n);
	for (int i = 1; i <= n; ++i) {
		ll x;
		cin >> x;
		sum.upd(i, i, x);
	}
	for (int i = 1; i <= m; ++i) {
		int op, x, y;
		cin >> op >> x >> y;
		if (op == 1) {
			ll val;
			cin >> val;
			sum.upd(x, y, val);
		}
		else {
			cout << sum.qry(x, y) << '\n';
		}
	}

	return 0;
}