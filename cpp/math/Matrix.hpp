#pragma once

#include "basics/Assert.hpp"
#include "utility/Identity.hpp"

template <class T, class Add = std::plus<>, class Mul = std::multiplies<>, auto Zero = nullptr>
struct Vector {
	int n;
	std::vector<T> a;

	Vector() : n() {}
	Vector(int n) : n(n), a(init(n)) {}
	Vector(int n, T v) : n(n), a(n, v) {}
	Vector(std::vector<T> v) : n(v.size()), a(std::move(v)) {}
	Vector(std::initializer_list<T> v) : n(v.size()), a(v) {}

	decltype(auto) operator[](int i) { return a[i]; }
	decltype(auto) operator[](int i) const { return a[i]; }

	Vector &operator+=(const Vector &b) {
		ASSERT(n == b.n);
		for (int i = 0; i < n; ++i) { add(a[i], b[i]); }
		return *this;
	}
	Vector &operator-=(const Vector &b) {
		ASSERT(n == b.n);
		for (int i = 0; i < n; ++i) { a[i] -= b[i]; }
		return *this;
	}
	Vector &operator*=(const T &k) {
		for (int i = 0; i < n; ++i) {
			if constexpr (std::is_same_v<Mul, std::multiplies<>> && !std::is_invocable_v<Mul, T&, const T&>) { a[i] *= k; }
			else { a[i] = Mul{}(a[i], k); }
		}
		return *this;
	}
	Vector &operator/=(const T &k) {
		for (int i = 0; i < n; ++i) { a[i] /= k; }
		return *this;
	}

	friend Vector operator+(Vector a, const Vector &b) { return a += b; }
	friend Vector operator-(Vector a, const Vector &b) { return a -= b; }
	friend Vector operator*(Vector a, const T &k) { return a *= k; }
	friend Vector operator*(const T &k, Vector a) {
		if constexpr (std::is_same_v<Mul, std::multiplies<>>) { a *= k; }
		else for (int i = 0; i < a.n; ++i) { a[i] = Mul{}(k, a[i]); }
		return a;
	}
	friend Vector operator/(Vector a, const T &k) { return a /= k; }

	T dot(const Vector &b) const {
		ASSERT(n == b.n);
		T res = zero();
		for (int i = 0; i < n; ++i) { add(res, Mul{}(a[i], b[i])); }
		return res;
	}

private:
	static constexpr T zero() { return identity_detail::get<T, Zero>(); }

	static std::vector<T> init(int n) {
		if constexpr (std::is_same_v<decltype(Zero), std::nullptr_t>) { return std::vector<T>(n); }
		else { return std::vector<T>(n, zero()); }
	}

	template <class L, class R>
	static void add(L &&x, R &&y) {
		if constexpr (std::is_same_v<Add, std::plus<>> && !std::is_invocable_v<Add, L&, R&&>) { x += std::forward<R>(y); }
		else { x = Add{}(x, std::forward<R>(y)); }
	}
};

template <class T, class Add = std::plus<>, class Mul = std::multiplies<>, auto Zero = nullptr, auto One = nullptr>
struct Matrix {
	int n, m;
	std::vector<std::vector<T>> a;

	Matrix() : n(), m() {}
	Matrix(int n, int m) : n(n), m(m), a(n, std::vector<T>(m, zero())) {}
	Matrix(int n, int m, T v) : Matrix(n, m) {
		for (int i = 0; i < n && i < m; ++i) { a[i][i] = v; }
	}
	Matrix(std::vector<std::vector<T>> v)
		: n(v.size()), m(n ? v[0].size() : 0), a(std::move(v)) {
		for (int i = 0; i < n; ++i) { ASSERT((int)a[i].size() == m); }
	}
	Matrix(std::initializer_list<std::initializer_list<T>> v)
		: n(v.size()), m(n ? v.begin()->size() : 0) {
		a.reserve(n);
		for (const auto &row : v) {
			ASSERT((int)row.size() == m);
			a.emplace_back(row);
		}
	}

	auto &operator[](int i) { return a[i]; }
	const auto &operator[](int i) const { return a[i]; }

	Matrix transpose() const {
		Matrix res(m, n);
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) { res[j][i] = a[i][j]; }
		}
		return res;
	}

	Matrix &operator+=(const Matrix &b) {
		ASSERT(n == b.n && m == b.m);
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) { a[i][j] = Add{}(a[i][j], b[i][j]); }
		}
		return *this;
	}

	Matrix &operator-=(const Matrix &b) {
		ASSERT(n == b.n && m == b.m);
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) { a[i][j] -= b[i][j]; }
		}
		return *this;
	}

	Matrix &operator*=(const Matrix &b) { return *this = *this * b; }
	friend Matrix operator+(Matrix a, const Matrix &b) { return a += b; }
	friend Matrix operator-(Matrix a, const Matrix &b) { return a -= b; }
	friend Matrix operator*(const Matrix &a, const Matrix &b) {
		ASSERT(a.m == b.n);
		Matrix res(a.n, b.m);
		for (int i = 0; i < a.n; ++i) {
			for (int k = 0; k < a.m; ++k) {
				for (int j = 0; j < b.m; ++j) {
					res[i][j] = Add{}(res[i][j], Mul{}(a[i][k], b[k][j]));
				}
			}
		}
		return res;
	}

	template <class A, class B, auto Z>
	friend auto operator*(const Matrix &a, const Vector<T, A, B, Z> &b) {
		ASSERT(a.m == b.n);
		Vector<T, A, B, Z> res(a.n, zero());
		for (int i = 0; i < a.n; ++i) {
			for (int j = 0; j < a.m; ++j) {
				res[i] = Add{}(res[i], Mul{}(a[i][j], b[j]));
			}
		}
		return res;
	}

	template <class A, class B, auto Z>
	friend auto operator*(const Vector<T, A, B, Z> &a, const Matrix &b) {
		ASSERT(a.n == b.n);
		Vector<T, A, B, Z> res(b.m, zero());
		for (int i = 0; i < b.n; ++i) {
			for (int j = 0; j < b.m; ++j) {
				res[j] = Add{}(res[j], Mul{}(a[i], b[i][j]));
			}
		}
		return res;
	}

	Matrix pow(long long k) const {
		ASSERT(n == m && k >= 0);
		Matrix res(n, n, one()), x = *this;
		while (k) {
			if (k & 1) { res *= x; }
			k >>= 1;
			if (k) { x *= x; }
		}
		return res;
	}

private:
	static constexpr T zero() { return identity_detail::get<T, Zero>(); }
	static constexpr T one() { return identity_detail::get<T, One, true>(); }
};
