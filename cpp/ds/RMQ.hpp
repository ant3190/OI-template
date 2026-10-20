#pragma once

#include "basics/Assert.hpp"

template <class T>
struct RMQ {
private:
	using Op = std::function<T(const T&, const T&)>;

public:
	RMQ() : n(-1) {}

	template <class It>
	RMQ(It s, It e, Op best = [](const T& x, const T& y) { return std::min(x, y); }) :
		best(std::move(best)) { build(s, e); }

	RMQ(const std::vector<T>& vec, Op best = [](const T& x, const T& y) { return std::min(x, y); }) :
		RMQ(vec.begin(), vec.end(), std::move(best)) { }

	T query(int l, int r) const {
		ASSERT(0 <= l && l <= r && r < n);
		int k = 31 - __builtin_clz(r - l + 1);
		return best(st[k][l], st[k][r - (1 << k) + 1]);
	}

private:
	int n;
	Op best;
	std::vector<std::vector<T>> st;

	template <class It>
	void build(It s, It e) {
		std::vector<T> a(s, e);
		n = a.size();
		ASSERT(n > 0);

		st.resize(32 - __builtin_clz(n));
		st[0] = std::move(a);

		for (int i = 1; i < (int)st.size(); ++i) {
			st[i].resize(n - (1 << i) + 1);
			for (int j = 0; j < (int)st[i].size(); ++j) {
				st[i][j] = best(
					st[i - 1][j],
					st[i - 1][j + (1 << (i - 1))]
				);
			}
		}
	}
};