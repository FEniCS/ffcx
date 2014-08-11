"""Tests of pure C++ expression formatting tools, not involving any UFL expressions."""

# TODO: Change methods of these classes to not take ufl expression as first argument

import uflacs
from uflacs.codeutils.cpp_expr_formatting_rules import CppFormattingRules

from ufl import as_ufl

def test_cpp_literal_formatter():
    fmt = CppFormattingRules()

    assert fmt.int_value(as_ufl(123)) == "123"
    assert fmt.float_value(as_ufl(1.0/3.0))[:7] == "0.33333"
    assert fmt.zero(as_ufl(0)) == "0.0"

    assert fmt.int_value(123) == "123"
    assert fmt.float_value(1.0/3.0)[:7] == "0.33333"
    assert fmt.zero(0) == "0.0"

    assert fmt.int_value("123") == "123"
    assert fmt.float_value(str(1.0/3.0))[:7] == "0.33333"
    assert fmt.zero("0") == "0.0"

def test_cpp_arithmetic_formatter():
    fmt = CppFormattingRules()

    assert fmt.sum(None, "1", "2") == "1 + 2"
    assert fmt.sum(None, "1", "2", "3") == "1 + 2 + 3"
    assert fmt.product(None, "1", "2") == "1 * 2"
    assert fmt.product(None, "1", "2", "3") == "1 * 2 * 3"
    assert fmt.division(None, "1", "2") == "1 / 2"

def test_cpp_cmath_formatter():
    fmt = CppFormattingRules()

    std = "std::"
    #std = ""

    assert fmt.power(None, "x", "y") == "{0}pow(x, y)".format(std)
    assert fmt.sqrt(None, "x") == "{0}sqrt(x)".format(std)
    assert fmt.ln(None, "x") == "{0}log(x)".format(std)
    assert fmt.exp(None, "x") == "{0}exp(x)".format(std)
    assert fmt.abs(None, "x") == "{0}abs(x)".format(std)
    #assert fmt.abs(None, "x") == "{0}fabs(x)".format(std)
    assert fmt.cos(None, "x") == "{0}cos(x)".format(std)
    assert fmt.sin(None, "x") == "{0}sin(x)".format(std)
    assert fmt.tan(None, "x") == "{0}tan(x)".format(std)
    assert fmt.cosh(None, "x") == "{0}cosh(x)".format(std)
    assert fmt.sinh(None, "x") == "{0}sinh(x)".format(std)
    assert fmt.tanh(None, "x") == "{0}tanh(x)".format(std)
    assert fmt.acos(None, "x") == "{0}acos(x)".format(std)
    assert fmt.asin(None, "x") == "{0}asin(x)".format(std)
    assert fmt.atan(None, "x") == "{0}atan(x)".format(std)
    #self.assertEqual(fmt.bessel_i(None, "n", "x"), "cyl_bessel_i(n, x)") # TODO
    #self.assertEqual(fmt.bessel_j(None, "n", "x"), "cyl_bessel_j(n, x)")
    #self.assertEqual(fmt.bessel_k(None, "n", "x"), "cyl_bessel_k(n, x)")
    #self.assertEqual(fmt.bessel_l(None, "n", "x"), "cyl_neumann(n, x)")

def test_cpp_conditional_formatter():
    fmt = CppFormattingRules()

    assert fmt.eq(None, "x", "y") == "x == y"
    assert fmt.ne(None, "x", "y") == "x != y"
    assert fmt.lt(None, "x", "y") == "x < y"
    assert fmt.gt(None, "x", "y") == "x > y"
    assert fmt.le(None, "x", "y") == "x <= y"
    assert fmt.ge(None, "x", "y") == "x >= y"
    assert fmt.and_condition(None, "x", "y") == "x && y"
    assert fmt.or_condition(None, "x", "y") == "x || y"
    assert fmt.not_condition(None, "x") == "!x"
