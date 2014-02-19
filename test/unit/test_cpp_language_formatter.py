"""Tests of pure C++ expression formatting tools, not involving any UFL expressions."""

# TODO: Change methods of these classes to not take ufl expression as first argument

def test_cpp_literal_formatter(self):
    fmt = CppLiteralFormatterRules()

    assert fmt.int_value(as_ufl(123)) == "123"
    assert fmt.float_value(as_ufl(1.0/3.0))[:7] == "0.33333"
    assert fmt.zero(as_ufl(0)) == "0"

    assert fmt.int_value(123) == "123"
    assert fmt.float_value(1.0/3.0)[:7] == "0.33333"
    assert fmt.zero(0) == "0"

    assert fmt.int_value("123") == "123"
    assert fmt.float_value(str(1.0/3.0))[:7] == "0.33333"
    assert fmt.zero("0") == "0"

def test_cpp_arithmetic_formatter(self):
    fmt = CppArithmeticFormatterRules()
    assert fmt.sum(None, "1", "2") == "1 + 2"
    assert fmt.sum(None, "1", "2", "3") == "1 + 2 + 3"
    assert fmt.product(None, "1", "2") == "1 * 2"
    assert fmt.product(None, "1", "2", "3") == "1 * 2 * 3"
    assert fmt.division(None, "1", "2") == "1 / 2"

def test_cpp_cmath_formatter(self):
    fmt = CppCmathFormatterRules()
    assert fmt.power(None, "x", "y") == "pow(x, y)"
    assert fmt.sqrt(None, "x") == "sqrt(x)"
    assert fmt.ln(None, "x") == "log(x)"
    assert fmt.exp(None, "x") == "exp(x)"
    assert fmt.abs(None, "x") == "fabs(x)"
    assert fmt.cos(None, "x") == "cos(x)"
    assert fmt.sin(None, "x") == "sin(x)"
    assert fmt.tan(None, "x") == "tan(x)"
    assert fmt.cosh(None, "x") == "cosh(x)"
    assert fmt.sinh(None, "x") == "sinh(x)"
    assert fmt.tanh(None, "x") == "tanh(x)"
    assert fmt.acos(None, "x") == "acos(x)"
    assert fmt.asin(None, "x") == "asin(x)"
    assert fmt.atan(None, "x") == "atan(x)"
    #self.assertEqual(fmt.bessel_i(None, "n", "x"), "cyl_bessel_i(n, x)") # TODO
    #self.assertEqual(fmt.bessel_j(None, "n", "x"), "cyl_bessel_j(n, x)")
    #self.assertEqual(fmt.bessel_k(None, "n", "x"), "cyl_bessel_k(n, x)")
    #self.assertEqual(fmt.bessel_l(None, "n", "x"), "cyl_neumann(n, x)")

def test_cpp_conditional_formatter(self):
    fmt = CppConditionalFormatterRules()
    assert fmt.eq(None, "x", "y") == "x == y"
    assert fmt.ne(None, "x", "y") == "x != y"
    assert fmt.lt(None, "x", "y") == "x < y"
    assert fmt.gt(None, "x", "y") == "x > y"
    assert fmt.le(None, "x", "y") == "x <= y"
    assert fmt.ge(None, "x", "y") == "x >= y"
    assert fmt.and_condition(None, "x", "y") == "x && y"
    assert fmt.or_condition(None, "x", "y") == "x || y"
    assert fmt.not_condition(None, "x") == "!x"
