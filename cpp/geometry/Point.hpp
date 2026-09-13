#pragma once

#include "basics/Assert.hpp"

template<class T> constexpr T eps_v = static_cast<T>(1e-9L);

template<class T>
int sgn(T x) { return (x > eps_v<T>) - (x < -eps_v<T>); }

template<class T> 
struct Point {
	using P = Point;

	T x, y;
	Point() : x(), y() {}
	Point(T x, T y) : x(x), y(y) {}

	P& operator+=(const P& p) { x += p.x, y += p.y; return *this; }
	P& operator-=(const P& p) { x -= p.x, y -= p.y; return *this; }
	bool operator<(const P& p) const { return std::tie(x, y) < std::tie(p.x, p.y); }
	bool operator==(const P& p) const { return std::tie(x, y) == std::tie(p.x, p.y); }
	bool operator!=(const P& p) const { return !(*this == p); }
	P operator+(const P& p) const { return P(x + p.x, y + p.y); }
	P operator-(const P& p) const { return P(x - p.x, y - p.y); }
	T operator^(const P& p) const { return x * p.y - y * p.x; }
	T operator*(const P& p) const { return x * p.x + y * p.y; }
	P operator*(T d) const { return P(x * d, y * d); }
	P operator/(T d) const { return P(x / d, y / d); }
	friend P operator*(T d, const P& p) { return p * d; }

	T len2() const { return x * x + y * y; }
	double len() const { return std::sqrt((double)len2()); }
	double angle() const { return std::atan2((double)y, (double)x); }
};

template<class T> T dist2(Point<T> a, Point<T> b) { return (b - a).len2(); }
template<class P> double dist(P a, P b) { return (b - a).len(); }

template<class P> void polar_sort(std::vector<P> &v, P o = P()) {
	auto half = [](const P& p) {
		return p.y < 0 || (p.y == 0 && p.x < 0);
	};
	std::sort(v.begin(), v.end(), [&](const P& a, const P& b) {
		P qa = a - o, qb = b - o;
		int ha = half(qa), hb = half(qb);
		if (ha != hb) { return ha < hb; }
		auto c = qa ^ qb;
		return c != 0 ? c > 0 : qa.len2() < qb.len2();
	});
}