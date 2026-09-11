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
