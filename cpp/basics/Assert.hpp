#pragma once

#ifdef LOCAL
#include <cassert>
#define ASSERT assert
#else
#define ASSERT(...) ((void)0)
#endif
