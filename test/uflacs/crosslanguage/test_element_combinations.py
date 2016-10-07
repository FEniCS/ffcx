#!/usr/bin/env py.test
# -*- coding: utf-8 -*-

class MockBasicElementCodeGenerator:
    pass

def test_mock_element_generator(gtest):
    """Test that the basic element mock class works as assumed."""

    pre = """
    double x[3] = { 0.0, 0.1, 0.2 };
    """

    post = """
    ASSERT_EQ(xi[0], 0.0);
    ASSERT_EQ(xi[1], 0.1);
    ASSERT_EQ(xi[2], 0.2);
    """

    # FIXME: Define interface for mock element code generator using TDD

    code = """
    double xi[3];
    xi[0] = x[0];
    xi[1] = x[1];
    xi[2] = x[2];
    """

    gtest.add(pre + code + post)

