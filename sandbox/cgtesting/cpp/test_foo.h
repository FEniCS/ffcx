#ifndef TEST_FOO_H_
#define TEST_FOO_H_

#include <gtest/gtest.h>

#include "gc_test_foo.h"

#include <vector>
#include <string>
using namespace std;

TEST (test_foo, test_some_generated_code)
{
    double v = some();
    ASSERT_EQ(v, 0.1);
}

TEST (test_foo, test_some_more_generated_code)
{
    double v = some_more();
    ASSERT_EQ(v, 0.2);
}

#endif
