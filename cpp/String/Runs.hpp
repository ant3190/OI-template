#pragma once

#include "String/SuffixArray.hpp"

template<class S>
std::vector<std::array<int, 3>> get_runs(const S& s) {
	int n = (int)s.size() - 1;
	if (n <= 1) { return {}; }
	SuffixArray sa(s);
	sa.init_lcp();
	S rs = s;
	std::reverse(rs.begin() + 1, rs.end());
	SuffixArray rsa(rs);
	rsa.init_lcp();

	std::vector<std::array<int, 3>> res;
	auto solve = [&](const std::vector<int>& rk) {
		std::vector<int> la(n + 1), stk;
		stk.reserve(n);
		for (int i = n; i >= 1; --i) {
			while (!stk.empty() && rk[stk.back()] > rk[i]) {
				stk.pop_back();
			}
			la[i] = stk.empty() ? n - i + 1 : stk.back() - i;
			stk.push_back(i);
		}
		for (int i = 1; i <= n; ++i) {
			int p = la[i];
			if (i + p > n) { continue; }
			int r = sa.lcp(i, i + p);
			int l = (i > 1 ? rsa.lcp(n - i + 2, n - i - p + 2) : 0);
			if (l + r >= p) {
				res.push_back({i - l, i + p + r - 1, p});
			}
		}
	};

	solve(sa.rk);
	std::vector<int> ord(n + 1);
	int m = (int)rank_compress(s.begin() + 1, s.end(), ord.begin() + 1).size();
	for (int i = 1; i <= n; ++i) { ord[i] = m + 1 - ord[i]; }
	SuffixArray isa(ord, m);
	solve(isa.rk);

	std::sort(res.begin(), res.end());
	res.erase(std::unique(res.begin(), res.end()), res.end());

	return res;
}