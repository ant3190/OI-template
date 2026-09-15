#pragma once

#include "basics/Assert.hpp"
#include "ds/RMQ.hpp"

template <class T, class Cmp = std::less<T>>
struct LinearRMQ {
public:
	LinearRMQ() : n(-1) {}
	template <class It>
	LinearRMQ(It s, It e) { build(s, e); }
	LinearRMQ(const std::vector<T> &vec) { build(vec.begin(), vec.end()); }

	T query(int l, int r) const {
		ASSERT(0 <= l && l <= r && r < n);
		if (r - l < B) { return small(l, r); }
		int bl = l >> LG, br = r >> LG;
		T res = best(small(l, (bl + 1) * B - 1), small(br * B, r));
		if (bl + 1 < br) {
			res = best(res, rmq.query(bl + 1, br - 1));
		}
		return res;
	}

private:
	static constexpr int B = 64, LG = 6;
	using ull = unsigned long long;

	int n;
	Cmp cmp;
	std::vector<T> a;
	std::vector<ull> mask;
	RMQ<T, Cmp> rmq;

	const T &best(const T &x, const T &y) const {
		return cmp(y, x) ? y : x;
	}

	T small(int l, int r) const {
		ull s = mask[r];
		int len = r - l + 1;
		if (len < B) { s &= (1ULL << len) - 1; }
		return a[r - (63 - __builtin_clzll(s))];
	}

	template <class It>
	void build(It s, It e) {
		a.assign(s, e);
		n = a.size();
		ASSERT(n > 0);

		mask.resize(n);
		for (int i = 0; i < n; ++i) {
			ull cur = (i ? mask[i - 1] << 1 : 0);
			while (cur) {
				int d = __builtin_ctzll(cur);
				if (cmp(a[i - d], a[i])) { break; }
				cur &= cur - 1;
			}
			mask[i] = cur | 1;
		}
		std::vector<T> b((n + B - 1) / B);
		for (int i = 0; i < (int)b.size(); ++i) {
			int l = i * B;
			int r = std::min(n, l + B) - 1;
			b[i] = small(l, r);
		}
		rmq = RMQ<T, Cmp>(b);
	}
};