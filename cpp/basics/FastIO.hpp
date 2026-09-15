#pragma once

namespace FastIO {
	static char buf[100000], *p1 = buf, *p2 = buf;
	#define gc (p1 == p2 && (p2 = (p1 = buf) + std::fread(buf, 1, 100000, stdin), p1 == p2) ? EOF : *p1++)
	inline long long read() { 
		long long res = 0;
		int w = 0, c = gc; 
		for (; !std::isdigit(c); c = gc) {
			((c == '-') && (w = 1));
		}
		for (; std::isdigit(c); c = gc) {
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
	inline std::string readS() {
		std::string res = "";
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
		for (; !std::isdigit(c); c = gc) {
			((c == '-') && (w = 1));
		}
		for (; std::isdigit(c); c = gc) {
			res = (res * 10) + (c ^ 48);
		}
		if (c == '.') {
			c = gc;
			for (; std::isdigit(c); c = gc) {
				res = res + tmp * (c ^ 48);
				tmp *= 0.1;
			}
		}
		return (w ? -res : res);
	}
	inline void write(long long x, char c = '\n') {
		((x < 0) && (std::putchar('-'), x *= -1));
		static int sta[50], top = 0; 
		do {
			sta[top++] = x % 10, x /= 10;
		} while (x); 
		while (top) {
			std::putchar(sta[--top] + 48);
		} 
		std::putchar(c);
	}
}

using namespace FastIO;