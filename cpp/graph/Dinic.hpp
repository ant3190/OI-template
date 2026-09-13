#pragma once

#include "basics/Assert.hpp"

template<class T>
struct Dinic {
public:
	Dinic() : n() {}
	Dinic(int n) : n(n), ed(n + 1), dis(n + 1), qu(n + 1), ptr(n + 1) {}

	std::pair<int, int> add(int from, int to, T cap) {
		assert(cap >= 0);
		int id = (int)ed[from].size(), rev = (int)ed[to].size() + (from == to);
		ed[from].push_back(edge(to, rev, cap));
		ed[to].push_back(edge(from, id, 0));
		return {from, id};
	}

	T max_flow(int s, int t, T lim = inf) {
		assert(s != t);
		T res = 0;
		while (res < lim && path(s, t)) {
			res += dfs(s, t, lim - res);
		}
		return res;
	}

	std::pair<T, T> get_flow(std::pair<int, int> id) {
		auto [u, i] = id;
		edge &e = ed[u][i], &r = ed[e.to][e.rev];
		return {r.cap, e.cap + r.cap};
	}

	bool left_of_cut(int s) {
		return dis[s] != -1;
	}

	void reserve(int u, int m) { ed[u].reserve(m); }

private:
	static constexpr T inf = std::numeric_limits<T>::max() / 2;

	struct edge {
		int to, rev;
		T cap;
		edge() {}
		edge(int to, int rev, T cap) : to(to), rev(rev), cap(cap) {}
	};

	int n;

	std::vector<std::vector<edge>> ed;
	std::vector<int> dis, qu, ptr;

	bool path(int s, int t) {
		std::fill(dis.begin(), dis.end(), -1);
		std::fill(ptr.begin(), ptr.end(), 0);
		dis[s] = 0;
		qu[0] = s;
		int hd = 0, tl = 1;
		while (hd < tl && dis[t] == -1) {
			int u = qu[hd++];
			for (edge& e : ed[u]) {
				if (e.cap && dis[e.to] == -1) {
					dis[e.to] = dis[u] + 1;
					qu[tl++] = e.to;
				}
			}
		}
		return dis[t] != -1;
	}

	T dfs(int u, int t, T flw) {
		if (u == t) { return flw; }
		T res = 0;
		for (int& i = ptr[u]; i < (int)ed[u].size(); ++i) {
			edge& e = ed[u][i];
			if (e.cap && dis[e.to] == dis[u] + 1) {
				T tmp = dfs(e.to, t, std::min(flw, e.cap));
				res += tmp, flw -= tmp;
				e.cap -= tmp, ed[e.to][e.rev].cap += tmp;
				if (!flw) { return res; }
			}
		}
		if (!res) { dis[u] = -1; }
		return res;
	}
};