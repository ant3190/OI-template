
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