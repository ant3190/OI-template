#pragma once

template<class It, class Ot>
auto rank_compress(It s, It e, Ot p) {
	using T = typename std::iterator_traits<It>::value_type;
	std::vector<T> v(s, e), a = v;
	std::sort(v.begin(), v.end());
	v.erase(std::unique(v.begin(), v.end()), v.end());
	for (auto x : a) { *p++ = std::lower_bound(v.begin(), v.end(), x) - v.begin() + 1; }
	return v;
}

template<class It>
auto rank_compress(It s, It e) {
	return rank_compress(s, e, s);
}