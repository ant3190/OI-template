#pragma once

#include "basics/Assert.hpp"
#include "utility/Identity.hpp"

template <class T, int N, class Add = std::plus<>, class Mul = std::multiplies<>, auto Zero = nullptr>
struct StaticVector {
	static_assert(N >= 0);
	static constexpr int n = N;
	using Storage = std::array<T, N>;
	Storage a = init();

	constexpr StaticVector() = default;
	explicit constexpr StaticVector(const T &v) : a(filled(v, std::make_index_sequence<N>{})) {}
	constexpr StaticVector(Storage v) : a(std::move(v)) {}
	constexpr StaticVector(std::initializer_list<T> v) {
		ASSERT((int)v.size() <= n);
		int i = 0;
		for (const auto &x : v) { a[i++] = x; }
	}

	constexpr T &operator[](int i) { return a[i]; }
	constexpr const T &operator[](int i) const { return a[i]; }

	constexpr StaticVector &operator+=(const StaticVector &b) {
		for (int i = 0; i < n; ++i) { add(a[i], b[i]); }
		return *this;
	}
	constexpr StaticVector &operator-=(const StaticVector &b) {
		for (int i = 0; i < n; ++i) { a[i] -= b[i]; }
		return *this;
	}
	constexpr StaticVector &operator*=(const T &k) {
		for (int i = 0; i < n; ++i) {
			if constexpr (std::is_same_v<Mul, std::multiplies<>> && !std::is_invocable_v<Mul, T&, const T&>) { a[i] *= k; }
			else { a[i] = Mul{}(a[i], k); }
		}
		return *this;
	}
	constexpr StaticVector &operator/=(const T &k) {
		for (int i = 0; i < n; ++i) { a[i] /= k; }
		return *this;
	}

	friend constexpr StaticVector operator+(StaticVector a, const StaticVector &b) { return a += b; }
	friend constexpr StaticVector operator-(StaticVector a, const StaticVector &b) { return a -= b; }
	friend constexpr StaticVector operator*(StaticVector a, const T &k) { return a *= k; }
	friend constexpr StaticVector operator*(const T &k, StaticVector a) {
		if constexpr (std::is_same_v<Mul, std::multiplies<>>) { a *= k; }
		else for (int i = 0; i < n; ++i) { a[i] = Mul{}(k, a[i]); }
		return a;
	}
	friend constexpr StaticVector operator/(StaticVector a, const T &k) { return a /= k; }

	constexpr T dot(const StaticVector &b) const {
		T res = zero();
		for (int i = 0; i < n; ++i) { add(res, Mul{}(a[i], b[i])); }
		return res;
	}

private:
	static constexpr T zero() { return identity_detail::get<T, Zero>(); }

	static constexpr Storage init() {
		if constexpr (std::is_same_v<decltype(Zero), std::nullptr_t>) { return {}; }
		else { return filled(zero(), std::make_index_sequence<N>{}); }
	}

	template <std::size_t... I>
	static constexpr Storage filled(const T &value, std::index_sequence<I...>) {
		return {{(static_cast<void>(I), value)...}};
	}

	template <class L, class R>
	static constexpr void add(L &&x, R &&y) {
		if constexpr (std::is_same_v<Add, std::plus<>> && !std::is_invocable_v<Add, L&, R&&>) { x += std::forward<R>(y); }
		else { x = Add{}(x, std::forward<R>(y)); }
	}
};

template <class T, int N, int M = N, class Add = std::plus<>, class Mul = std::multiplies<>, auto Zero = nullptr, auto One = nullptr>
struct StaticMatrix {
	static_assert(N >= 0 && M >= 0);
	static constexpr int n = N, m = M;
	using Row = std::array<T, M>;
	using Storage = std::array<Row, N>;
	Storage a = init();

	constexpr StaticMatrix() = default;
	explicit constexpr StaticMatrix(T v) : StaticMatrix() {
		for (int i = 0; i < n && i < m; ++i) { a[i][i] = v; }
	}
	constexpr StaticMatrix(Storage v) : a(std::move(v)) {}
	constexpr StaticMatrix(std::initializer_list<std::initializer_list<T>> v) {
		ASSERT((int)v.size() == n);
		int i = 0;
		for (const auto &row : v) {
			ASSERT((int)row.size() == m);
			int j = 0;
			for (const auto &x : row) { a[i][j++] = x; }
			++i;
		}
	}

	constexpr auto &operator[](int i) { return a[i]; }
	constexpr const auto &operator[](int i) const { return a[i]; }

	constexpr auto transpose() const {
		StaticMatrix<T, M, N, Add, Mul, Zero, One> res;
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) { res[j][i] = a[i][j]; }
		}
		return res;
	}

	constexpr StaticMatrix &operator+=(const StaticMatrix &b) {
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) { a[i][j] = Add{}(a[i][j], b[i][j]); }
		}
		return *this;
	}

	constexpr StaticMatrix &operator-=(const StaticMatrix &b) {
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) { a[i][j] -= b[i][j]; }
		}
		return *this;
	}

	constexpr StaticMatrix &operator*=(const StaticMatrix<T, M, M, Add, Mul, Zero, One> &b) { return *this = *this * b; }

	friend constexpr StaticMatrix operator+(StaticMatrix a, const StaticMatrix &b) { return a += b; }
	friend constexpr StaticMatrix operator-(StaticMatrix a, const StaticMatrix &b) { return a -= b; }
	template <int K>
	friend constexpr auto operator*(const StaticMatrix &a, const StaticMatrix<T, M, K, Add, Mul, Zero, One> &b) {
		StaticMatrix<T, N, K, Add, Mul, Zero, One> res;
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
	friend constexpr auto operator*(const StaticMatrix &a, const StaticVector<T, M, A, B, Z> &b) {
		StaticVector<T, N, A, B, Z> res(zero());
		for (int i = 0; i < a.n; ++i) {
			for (int j = 0; j < a.m; ++j) {
				res[i] = Add{}(res[i], Mul{}(a[i][j], b[j]));
			}
		}
		return res;
	}

	template <class A, class B, auto Z>
	friend constexpr auto operator*(const StaticVector<T, N, A, B, Z> &a, const StaticMatrix &b) {
		StaticVector<T, M, A, B, Z> res(zero());
		for (int i = 0; i < b.n; ++i) {
			for (int j = 0; j < b.m; ++j) {
				res[j] = Add{}(res[j], Mul{}(a[i], b[i][j]));
			}
		}
		return res;
	}

	constexpr StaticMatrix pow(long long k) const {
		static_assert(n == m);
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

	template <class U, std::size_t... I>
	static constexpr auto fill(const U &v, std::index_sequence<I...>) {
		return std::array<U, sizeof...(I)>{{(static_cast<void>(I), v)...}};
	}
	static constexpr Storage init() {
		const Row row = fill(zero(), std::make_index_sequence<M>{});
		return fill(row, std::make_index_sequence<N>{});
	}
};
