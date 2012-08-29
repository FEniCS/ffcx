/*******************************
 *
 *  Main file for running all C++ tests in one go
 *
 ********************************/

#include <gtest/gtest.h>

#include "main_testall.h"

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
