#pragma once

#include "basics/Assert.hpp"

template<class S>
std::vector<std::pair<int, int>> duval(const S& s) {
	int n = (int)s.size() - 1;

	std::vector<std::pair<int, int>> res;
	for (int i = 1; i <= n;) {
		int j = i, k = i + 1;
		for (; k <= n && s[j] <= s[k]; ++k) {
			if (s[j] == s[k]) {
				++j;
			} else {
				j = i;
			}
		}
		while (i <= j) {
			res.emplace_back(i, i + k - j - 1);
			i += k - j;
		}
	}
	return res;
}