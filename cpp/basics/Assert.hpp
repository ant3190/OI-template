#pragma once

#ifdef LOCAL
#define ASSERT assert
#else
#define ASSERT(...) ((void)0)
#endif
