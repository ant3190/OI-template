#pragma once

#ifdef LOCAL
#define ASSERT(...) assert((__VA_ARGS__))
#else
#define ASSERT(...) ((void)0)
#endif

template <class T>
struct MCMF {
public:
	MCMF() : n(), dirty() {}
	MCMF(int n) : n((ASSERT(n >= 0), n)), dirty(), ed(n + 1), par(n + 1), dis(n + 1), pi(n + 1) {}

	std::pair<int, int> add(int from, int to, T cap, T cost) {
		ASSERT(1 <= from && from <= n);
		ASSERT(1 <= to && to <= n);
		ASSERT(cap >= 0);
		dirty |= (cap && pi[from] != inf && pi[from] + cost < pi[to]);
		int id = (int)ed[from].size(), rev = (int)ed[to].size() + (from == to);
		ed[from].emplace_back(to, rev, cap, cost);
		ed[to].emplace_back(from, id, 0, -cost);
		return {from, id};
	}

	template <class F>
	std::pair<int, int> add_convex(int from, int to, T cap, F cost) {
		ASSERT(cap >= 0);
		T val = cost(T(0));
		auto id = add(from, to, cap, cap ? cost(T(1)) - val : T(0));
		fun.push_back({std::function<T(T)>(std::move(cost)), val});
		auto [u, i] = id;
		ed[u][i].fn = (int)fun.size();
		ed[ed[u][i].to][ed[u][i].rev].fn = -(int)fun.size();
		return id;
	}

	std::pair<T, T> max_flow(int s, int t, T lim = inf) {
		ASSERT(1 <= s && s <= n);
		ASSERT(1 <= t && t <= n);
		ASSERT(s != t && lim >= 0);
		setpi(s);
		std::pair<T, T> res = {0, 0};
		while (res.first < lim && path(s, t)) {
			T flw = lim - res.first;
			for (int i = t, j, r; i != s; i = j) {
				r = ed[i][par[i]].rev, j = ed[i][par[i]].to;
				edge &e = ed[j][r];
				flw = std::min(flw, e.cap);
				if (e.fn || flw == 1) { flw = 1; break; }
			}
			res.first += flw, res.second += flw * (pi[t] - pi[s]);
			for (int i = t, j, r; i != s; i = j) {
				r = ed[i][par[i]].rev, j = ed[i][par[i]].to;
				ed[i][par[i]].cap += flw;
				ed[j][r].cap -= flw;
				if (ed[j][r].fn) { update(j, r); }
			}
		}
		return res;
	}

	std::pair<T, T> get_flow(std::pair<int, int> id) {
		auto [u, i] = id;
		ASSERT(1 <= u && u <= n);
		ASSERT(0 <= i && i < (int)ed[u].size());
		edge &e = ed[u][i], &r = ed[e.to][e.rev];
		return {r.cap, e.cap + r.cap};
	}

	void reserve(int u, int m) {
		ASSERT(1 <= u && u <= n && m >= 0);
		ed[u].reserve(m);
	}

private:
	static constexpr T inf = std::numeric_limits<T>::max() / 2;

	struct edge {
		int to, rev, fn;
		T cap, cost;
		edge(int to, int rev, T cap, T cost) : to(to), rev(rev), fn(), cap(cap), cost(cost) {}
	};

	struct convex {
		std::function<T(T)> cost;
		T sum;
	};

	int n;
	bool dirty;

	std::vector<std::vector<edge>> ed;
	std::vector<int> par;
	std::vector<T> dis, pi;
	std::vector<convex> fun;
	std::vector<std::pair<T, int>> heap;

	void update(int u, int i) {
		edge &e = ed[u][i], &r = ed[e.to][e.rev];
		convex &c = fun[(e.fn > 0 ? e.fn : -e.fn) - 1];
		T old = e.cost;
		c.sum += old;
		r.cost = -old;
		if (e.cap) {
			T x = e.fn > 0 ? r.cap + 1 : e.cap - 1;
			e.cost = c.cost(x) - c.sum;
		}
	}

	bool path(int s, int t) {
		std::fill(dis.begin(), dis.end(), inf);
		dis[s] = 0;
		heap.clear();
		heap.push_back({T(0), s});
		const auto cmp = std::greater<std::pair<T, int>>();
		while (!heap.empty()) {
			std::pop_heap(heap.begin(), heap.end(), cmp);
			auto [di, u] = heap.back(); heap.pop_back();
			if (di != dis[u]) { continue; }
			if (u == t) { break; }
			for (edge &e : ed[u]) {
				if (!e.cap) { continue; }
				T val = di + e.cost + pi[u] - pi[e.to];
				if (dis[e.to] > val) {
					dis[e.to] = val;
					par[e.to] = e.rev;
					heap.push_back({val, e.to});
					std::push_heap(heap.begin(), heap.end(), cmp);
				}
			}
		}
		if (dis[t] == inf) { return false; }
		for (int i = 1; i <= n; ++i) {
			if (pi[i] != inf) { pi[i] += std::min(dis[i], dis[t]); }
		}
		return true;
	}

	void setpi(int s) {
		if (!dirty && pi[s] != inf) { return ; }
		dirty = 0;
		std::fill(pi.begin(), pi.end(), inf);
		pi[s] = 0;
		int t = n, flg = 1;
		while (flg-- && t--) {
			for (int i = 1; i <= n; ++i) {
				if (pi[i] != inf) {
					for (edge &e : ed[i]) {
						if (e.cap && pi[e.to] > pi[i] + e.cost) {
							pi[e.to] = pi[i] + e.cost, flg = 1;
						}
					}
				}
			}
		}
		ASSERT(t >= 0);
	}
};
