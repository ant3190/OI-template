#pragma once

#include "basics/Assert.hpp"

template <class T>
struct Matrix {
	int n, m;
	std::vector<std::vector<T>> a;

	Matrix() : n(), m() {}
	Matrix(int n, int m) : n(n), m(m), a(n, std::vector<T>(m)) {}
	Matrix(int n, int m, T v) : Matrix(n, m) {
		for (int i = 0; i < std::min(n, m); ++i) {
			a[i][i] = v;
		}
	}
	Matrix(std::vector<std::vector<T>> v)
		: n(v.size()), m(n ? v[0].size() : 0), a(std::move(v)) {
		for (int i = 0; i < n; ++i) {
			assert((int)a[i].size() == m);
		}
	}

	std::vector<T> &operator[](int i) { return a[i]; }
	const std::vector<T> &operator[](int i) const { return a[i]; }

	Matrix transpose() const {
		Matrix res(m, n);
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				res[j][i] = a[i][j];
			}
		}
		return res;
	}

	Matrix &operator+=(const Matrix &b) {
		assert(n == b.n && m == b.m);
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				a[i][j] += b[i][j];
			}
		}
		return *this;
	}

	Matrix &operator-=(const Matrix &b) {
		assert(n == b.n && m == b.m);
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				a[i][j] -= b[i][j];
			}
		}
		return *this;
	}

	Matrix &operator*=(const Matrix &b) {
		assert(m == b.n);
		Matrix res(n, b.m);
		for (int i = 0; i < n; ++i) {
			for (int k = 0; k < m; ++k) {
				for (int j = 0; j < b.m; ++j) {
					res[i][j] += a[i][k] * b[k][j];
				}
			}
		}
		return *this = std::move(res);
	}

	friend Matrix operator+(Matrix a, const Matrix &b) { return a += b; }
	friend Matrix operator-(Matrix a, const Matrix &b) { return a -= b; }
	friend Matrix operator*(Matrix a, const Matrix &b) { return a *= b; }

	Matrix pow(long long k) const {
		assert(n == m && k >= 0);
		Matrix res(n, n, T(1)), x = *this;
		for (; k; k >>= 1, x *= x) {
			if (k & 1) {
				res *= x;
			}
		}
		return res;
	}
};
