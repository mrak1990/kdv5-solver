#ifndef KDV5_GTEST_EXT_H
#define KDV5_GTEST_EXT_H

#include "gtest/gtest.h"

namespace testing
{
    namespace internal
    {
        enum GTestColor {
            COLOR_DEFAULT,
            COLOR_RED,
            COLOR_GREEN,
            COLOR_YELLOW
        };

        extern void ColoredPrintf(GTestColor color, const char* fmt, ...);
    }
}
#define PRINTF(...) ((void)0)
#define _PRINTF(...) do { testing::internal::ColoredPrintf(testing::internal::COLOR_GREEN, "[          ] "); testing::internal::ColoredPrintf(testing::internal::COLOR_YELLOW, __VA_ARGS__); } while(0)

#endif //KDV5_GTEST_EXT_H
