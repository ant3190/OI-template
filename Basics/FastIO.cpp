#include <bits/stdc++.h>
#include <cassert>
#define INF 0x3f3f3f3f
#define INF_LL 0x3f3f3f3f3f3f3f3f
using namespace std;

using ll = long long;
using ull = unsigned long long;
using i128 = __int128_t;

mt19937 rd(chrono::steady_clock().now().time_since_epoch().count());
mt19937_64 rdll(chrono::steady_clock().now().time_since_epoch().count());

/*========== Begin ==========*/

namespace FastIO {
	static char buf[100000], *p1 = buf, *p2 = buf;
	#define gc (p1 == p2 && (p2 = (p1 = buf) + fread(buf, 1, 100000, stdin), p1 == p2) ? EOF : *p1++)
	inline ll read() { 
		ll res = 0;
		int w = 0, c = gc; 
		for (; !isdigit(c); c = gc) {
			((c == '-') && (w = 1));
		}
		for (; isdigit(c); c = gc) {
			res = (res << 1) + (res << 3) + (c ^ 48);
		}
		return (w ? -res : res);
	}
	inline char readC() {  
		int c = gc; 
		while (c == '\n' || c == '\r' || c == ' ') {
			c = gc;
		}
		return c;
	}
	inline string readS() {
		string res = "";
		char c = gc; 
		for (; (c == '\n' || c == '\r' || c == ' ' || c == EOF); c = gc);
		for (; !(c == '\n' || c == '\r' || c == ' ' || c == EOF); c = gc) {
			res += c;
		}
		return res;
	}	
	inline double readF() { 
		double res = 0, tmp = 0.1;
		int w = 0; 
		char c = gc; 
		for (; !isdigit(c); c = gc) {
			((c == '-') && (w = 1));
		}
		for (; isdigit(c); c = gc) {
			res = (res * 10) + (c ^ 48);
		}
		if (c == '.') {
			c = gc;
			for (; isdigit(c); c = gc) {
				res = res + tmp * (c ^ 48);
				tmp *= 0.1;
			}
		}
		return (w ? -res : res);
	}
	inline void write(ll x, char c = '\n') {
		((x < 0) && (putchar('-'), x *= -1));
		static int sta[50], top = 0; 
		do {
			sta[top++] = x % 10, x /= 10;
		} while (x); 
		while (top) {
			putchar(sta[--top] + 48);
		} 
		putchar(c);
	}
}

using namespace FastIO;

/*========== End ==========*/

int main() {
#ifdef LOCAL
	assert(freopen("test.in", "r", stdin));
	assert(freopen("test.out", "w", stdout));
#endif

	int n = read();
	ull res = 0;
	for (int i = 1; i <= n; ++i) {
		res ^= read();
	}
	write(res);

	return 0;
}