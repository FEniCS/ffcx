#ifndef TEST_BAR_H_INCLUDED
#define TEST_BAR_H_INCLUDED

#include <gtest/gtest.h>

#include "gc_test_bar.h"

TEST (test_bar, test_something)
{
    double v = some2() + some_more2();
    ASSERT_EQ(v, 0.1+0.2);
}

#endif
