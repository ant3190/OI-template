#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <utility>
#include "basics/Assert.hpp"
#include "utility/Identity.hpp"

namespace static_matrix_detail {

template <class T, std::size_t... I>
constexpr std::array<T, sizeof...(I)> repeated(const T &value, std::index_sequence<I...>) {
	return {{(static_cast<void>(I), value)...}};
}

template <std::size_t N, class T>
constexpr std::array<T, N> filled(const T &value) {
	return repeated(value, std::make_index_sequence<N>{});
}

}

template <class T, int N>
struct StaticVector {
	static_assert(N >= 0);
	static constexpr int n = N;
	using Storage = std::array<T, N>;
	Storage a{};

	constexpr StaticVector() = default;
	explicit constexpr StaticVector(const T &v) : a(static_matrix_detail::filled<N>(v)) {}
	explicit constexpr StaticVector(Storage v) : a(std::move(v)) {}
	constexpr StaticVector(std::initializer_list<T> v) {
		ASSERT(v.size() <= static_cast<std::size_t>(N));
		int i = 0;
		for (const auto &x : v) { a[i++] = x; }
	}

	constexpr T &operator[](int i) { return a[i]; }
	constexpr const T &operator[](int i) const { return a[i]; }

	constexpr StaticVector &operator+=(const StaticVector &b) {
		for (int i = 0; i < N; ++i) { a[i] += b[i]; }
		return *this;
	}
	constexpr StaticVector &operator-=(const StaticVector &b) {
		for (int i = 0; i < N; ++i) { a[i] -= b[i]; }
		return *this;
	}
	constexpr StaticVector &operator*=(const T &k) {
		for (auto &x : a) { x *= k; }
		return *this;
	}
	constexpr StaticVector &operator/=(const T &k) {
		for (auto &x : a) { x /= k; }
		return *this;
	}

	friend constexpr StaticVector operator+(StaticVector a, const StaticVector &b) { return a += b; }
	friend constexpr StaticVector operator-(StaticVector a, const StaticVector &b) { return a -= b; }
	friend constexpr StaticVector operator*(StaticVector a, const T &k) { return a *= k; }
	friend constexpr StaticVector operator*(const T &k, StaticVector a) { return a *= k; }
	friend constexpr StaticVector operator/(StaticVector a, const T &k) { return a /= k; }

	constexpr T dot(const StaticVector &b) const {
		T res = T();
		for (int i = 0; i < N; ++i) { res += a[i] * b[i]; }
		return res;
	}
};

template <class T, int R, int C = R, class Add = std::plus<>, class Mul = std::multiplies<>,
		auto Zero = nullptr, auto One = nullptr>
struct StaticMatrix {
	static_assert(R >= 0 && C >= 0);
	static constexpr int n = R, m = C;
	using Row = std::array<T, C>;
	using Storage = std::array<Row, R>;
	Storage a = static_matrix_detail::filled<R>(static_matrix_detail::filled<C>(zero()));

	constexpr StaticMatrix() = default;
	explicit constexpr StaticMatrix(T v) {
		for (int i = 0; i < std::min(R, C); ++i) {
			a[i][i] = v;
		}
	}
	explicit constexpr StaticMatrix(Storage v) : a(std::move(v)) {}

	constexpr Row &operator[](int i) { return a[i]; }
	constexpr const Row &operator[](int i) const { return a[i]; }

	constexpr StaticMatrix<T, C, R, Add, Mul, Zero, One> transpose() const {
		StaticMatrix<T, C, R, Add, Mul, Zero, One> res;
		for (int i = 0; i < R; ++i) {
			for (int j = 0; j < C; ++j) {
				res[j][i] = a[i][j];
			}
		}
		return res;
	}

	constexpr StaticMatrix &operator+=(const StaticMatrix &b) {
		for (int i = 0; i < R; ++i) {
			for (int j = 0; j < C; ++j) {
				a[i][j] = Add{}(a[i][j], b[i][j]);
			}
		}
		return *this;
	}

	constexpr StaticMatrix &operator-=(const StaticMatrix &b) {
		for (int i = 0; i < R; ++i) {
			for (int j = 0; j < C; ++j) {
				a[i][j] -= b[i][j];
			}
		}
		return *this;
	}

	constexpr StaticMatrix &operator*=(const StaticMatrix<T, C, C, Add, Mul, Zero, One> &b) {
		return *this = *this * b;
	}

	friend constexpr StaticMatrix operator+(StaticMatrix a, const StaticMatrix &b) { return a += b; }
	friend constexpr StaticMatrix operator-(StaticMatrix a, const StaticMatrix &b) { return a -= b; }
	template <int K>
	friend constexpr StaticMatrix<T, R, K, Add, Mul, Zero, One> operator*(
			const StaticMatrix &a, const StaticMatrix<T, C, K, Add, Mul, Zero, One> &b) {
		StaticMatrix<T, R, K, Add, Mul, Zero, One> res;
		for (int i = 0; i < R; ++i) {
			for (int k = 0; k < C; ++k) {
				for (int j = 0; j < K; ++j) {
					res[i][j] = Add{}(res[i][j], Mul{}(a[i][k], b[k][j]));
				}
			}
		}
		return res;
	}

	friend constexpr StaticVector<T, R> operator*(const StaticMatrix &a, const StaticVector<T, C> &b) {
		StaticVector<T, R> res(zero());
		for (int i = 0; i < R; ++i) {
			for (int j = 0; j < C; ++j) {
				res[i] = Add{}(res[i], Mul{}(a[i][j], b[j]));
			}
		}
		return res;
	}

	friend constexpr StaticVector<T, C> operator*(const StaticVector<T, R> &a, const StaticMatrix &b) {
		StaticVector<T, C> res(zero());
		for (int i = 0; i < R; ++i) {
			for (int j = 0; j < C; ++j) {
				res[j] = Add{}(res[j], Mul{}(a[i], b[i][j]));
			}
		}
		return res;
	}

	constexpr StaticMatrix pow(long long k) const {
		static_assert(R == C);
		ASSERT(k >= 0);
		StaticMatrix res(one()), x = *this;
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
