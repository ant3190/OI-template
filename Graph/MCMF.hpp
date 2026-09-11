template<class T>
struct MCMF {
public:
	MCMF() : n(), dirty() {}
	MCMF(int n) : n(n), ed(n + 1), par(n + 1), dis(n + 1), pi(n + 1), dirty() {}

	pair<int, int> adde(int from, int to, T cap, T cost) {
		assert(cap >= 0);
		dirty |= (cap && pi[from] != inf && pi[from] + cost < pi[to]);
		int id = (int)ed[from].size(), rev = (int)ed[to].size() + (from == to);
		ed[from].push_back(edge(to, rev, cap, cost));
		ed[to].push_back(edge(from, id, 0, -cost));
		return {from, id};
	}

	pair<T, T> max_flow(int s, int t, T lim = inf) {
		assert(s != t);
		setpi(s);
		pair<T, T> res = {0, 0};
		while (res.first < lim && path(s, t)) {
			T flw = lim - res.first;
			for (int i = t, j, r; i != s; i = j) {
				r = ed[i][par[i]].rev, j = ed[i][par[i]].to;
				flw = min(flw, ed[j][r].cap);
			}
			res.first += flw, res.second += flw * (pi[t] - pi[s]);
			for (int i = t, j, r; i != s; i = j) {
				r = ed[i][par[i]].rev, j = ed[i][par[i]].to;
				ed[i][par[i]].cap += flw;
				ed[j][r].cap -= flw;
			}
		}
		return res;
	}

	pair<T, T> getflow(pair<int, int> id) {
		auto [u, i] = id;
		edge &e = ed[u][i], &r = ed[e.to][e.rev];
		return {r.cap, e.cap + r.cap};
	}

	void reserve(int u, int m) { ed[u].reserve(m); }

private:
	static constexpr T inf = numeric_limits<T>::max() / 2;

	struct edge {
		int to, rev;
		T cap, cost;
		edge() {}
		edge(int to, int rev, T cap, T cost) : to(to), rev(rev), cap(cap), cost(cost) {}
	};

	int n;
	bool dirty;

	vector<vector<edge>> ed;
	vector<int> par;
	vector<T> dis, pi;

	bool path(int s, int t) {
		fill(dis.begin(), dis.end(), inf);
		dis[s] = 0;
		priority_queue<pair<T, int>, vector<pair<T, int>>, greater<pair<T, int>>> qu;
		qu.push({dis[s], s});
		while (!qu.empty()) {
			auto [di, u] = qu.top(); qu.pop();
			if (di != dis[u]) { continue; }
			for (edge& e : ed[u]) {
				T val = di + e.cost + pi[u] - pi[e.to];
				if (e.cap && dis[e.to] > val) {
					dis[e.to] = val;
					par[e.to] = e.rev;
					qu.push({dis[e.to], e.to});
				}
			}
		}
		for (int i = 1; i <= n; ++i) {
			pi[i] = (dis[i] == inf ? inf : pi[i] + dis[i]);
		}
		return dis[t] != inf;
	}

	void setpi(int s) {
		if (!dirty && pi[s] != inf) { return ; }
		dirty = 0;

		fill(pi.begin(), pi.end(), inf);
		pi[s] = 0;
		int t = n, flg = 1;
		while (flg-- && t--) {
			for (int i = 1; i <= n; ++i) {
				if (pi[i] != inf) {
					for (auto& e : ed[i]) {
						if (e.cap && pi[e.to] > pi[i] + e.cost) {
							pi[e.to] = pi[i] + e.cost, flg = 1;
						}
					}
				}
			}
		}
		assert(t >= 0);
	}
};