#!/usr/bin/env py.test
# -*- coding: utf-8 -*-

import ufl
from ufl import product

from ffc.uflacs.backends.toy.toy_compiler import compile_expression

"""
Unit tests of generated geometry snippet code.
"""

def compile_expression0(expr):
    code = ""
    return code

def compile_expression1(expr):
    code = "double values[%d];" % (product(expr.ufl_shape),)
    return code

def compile_expression2(expr):
    code = compile_expression(expr, "")
    return code

def test_compilation_of_x(gtest):
    """Test that compilation of x results in the value of x placed in values."""

    pre = """
    double x[3] = { 0.0, 0.1, 0.2 };
    double A[3];
    """

    post = """
    ASSERT_EQ(A[0], 0.0);
    ASSERT_EQ(A[1], 0.1);
    ASSERT_EQ(A[2], 0.2);
    """

    expr = ufl.SpatialCoordinate(ufl.tetrahedron)
    code = compile_expression2(expr)

    gtest.add(pre + code + post)

def test_compilation_of_sums(gtest):
    """Test that sums are compiled correctly."""

    pre = """
    double x[3] = { 0.1, 0.2, 0.3 };
    double A[1];
    """

    post = """
    ASSERT_EQ(A[0], x[0] + x[1] + x[2]);
    """

    x = ufl.SpatialCoordinate(ufl.tetrahedron)
    expr = x[0] + x[1] + x[2]
    code = compile_expression2(expr)
    print(code)

    gtest.add(pre + code + post)

def test_compilation_of_products(gtest):
    """Test that products are compiled correctly."""

    pre = """
    double x[3] = { 0.1, 0.2, 0.3 };
    double A[1];
    """

    post = """
    ASSERT_EQ(A[0], x[0] * x[1] * x[2]);
    """

    x = ufl.SpatialCoordinate(ufl.tetrahedron)
    expr = x[0] * x[1] * x[2]
    code = compile_expression2(expr)
    print(code)

    gtest.add(pre + code + post)

def test_compilation_of_sums_and_products_with_precedence(gtest):
    """Test that combinations of sums and products are
    compiled correctly with precedence rules in mind.
    """

    pre = """
    double x[3] = { 0.1, 0.2, 0.3 };
    double A[1];
    """

    post = """
    ASSERT_EQ(A[0], (x[0] + x[1]) * x[2] - (x[0] + (x[1] * x[2])));
    """

    x = ufl.SpatialCoordinate(ufl.tetrahedron)
    expr = (x[0] + x[1]) * x[2] - (x[0] + (x[1] * x[2]))
    code = compile_expression2(expr)
    print(code)

    gtest.add(pre + code + post)

