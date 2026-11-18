#pragma once

template <class Cmp = std::less<>>
struct Best {
    template <class T>
    constexpr T operator()(const T& x, const T& y) const {
        return Cmp{}(y, x) ? y : x;
    }
};