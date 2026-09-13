#pragma once

#ifdef LOCAL
#define ASSERT(...) assert((__VA_ARGS__))
#else
#define ASSERT(...) ((void)0)
#endif
