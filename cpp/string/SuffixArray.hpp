#pragma once

#include "basics/Assert.hpp"
#include "utility/Order.hpp"

struct SuffixArray {
public:
	std::vector<int> sa, rk, ht;

	SuffixArray() : n(-1) {};
	template<class S> 
	SuffixArray(const S& s, int m = -1) {
		ASSERT(s.size() > 1);
		n = (int)s.size() - 1;
		sa.resize(n + 1), rk.resize(n + 1), ht.resize(n + 1);

		if (m == -1) {
			m = rank_compress(s.begin() + 1, s.end(), rk.begin() + 1).size();
		} else {
			for (int i = 1; i <= n; ++i) {
				rk[i] = s[i];
				ASSERT(1 <= rk[i] && rk[i] <= m);
			}
		}

		std::vector<int> cnt(std::max(m, n) + 1, 0), aux(n + 1), nrk(n + 1);
		for (int i = 1; i <= n; ++i) { ++cnt[rk[i]]; }
		for (int i = 1; i <= m; ++i) { cnt[i] += cnt[i - 1]; }
		for (int i = n; i >= 1; --i) { sa[cnt[rk[i]]--] = i; }
		for (int k = 1; k == 1 || m != n; k <<= 1) {
			int p = 0;
			for (int i = n - k + 1; i <= n; ++i) { aux[++p] = i; }
			for (int i = 1; i <= n; ++i) {
				if (sa[i] > k) { aux[++p] = sa[i] - k; }
			}
			std::fill(cnt.begin() + 1, cnt.begin() + m + 1, 0);
			for (int i = 1; i <= n; ++i) { ++cnt[rk[i]]; }
			for (int i = 1; i <= m; ++i) { cnt[i] += cnt[i - 1]; }
			for (int i = n; i >= 1; --i) { sa[cnt[rk[aux[i]]]--] = aux[i]; }
			nrk[sa[1]] = m = 1;
			for (int i = 2; i <= n; ++i) {
				int a = sa[i], b = sa[i - 1];
				m += (rk[a] != rk[b] || (a + k <= n ? rk[a + k] : -1) != (b + k <= n ? rk[b + k] : -1));
				nrk[a] = m;
			}
			std::swap(nrk, rk);
		}

		for (int i = 1, j, k = 0; i <= n; ++i) {
			if (rk[i] == 1) { continue; }
			if (k) { --k; }
			j = sa[rk[i] - 1];
			while (i + k <= n && j + k <= n && s[i + k] == s[j + k]) { ++k; }
			ht[rk[i]] = k;
		}
	}

	int size() const { return n; }
	int operator[](const int &id) const {
		ASSERT(id > 0 && id <= n);
		return sa[id];
	}

	void init_lcp() {
		ASSERT(n != -1);
		int m = 32 - __builtin_clz(n);
		rmq = std::vector<std::vector<int>>(m, std::vector<int>(n + 1));
		for (int i = 2; i <= n; ++i) { rmq[0][i] = ht[i]; }
		for (int i = 1; i < m; ++i) {
			for (int j = 2; j + (1 << i) - 1 <= n; ++j) {
				rmq[i][j] = std::min(rmq[i - 1][j], rmq[i - 1][j + (1 << (i - 1))]);
			}
		}
	}

	int lcp(int i, int j) {
		ASSERT(!rmq.empty());
		ASSERT(1 <= i && i <= n);
		ASSERT(1 <= j && j <= n);
		if (i == j) { return n - i + 1; }
		i = rk[i], j = rk[j];
		if (i > j) { std::swap(i, j); }
		int t = 31 - __builtin_clz(j - i++);
		return std::min(rmq[t][i], rmq[t][j - (1 << t) + 1]);
	}

private:
	int n;
	std::vector<std::vector<int>> rmq;
};