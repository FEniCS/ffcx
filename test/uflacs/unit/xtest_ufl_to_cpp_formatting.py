# -*- coding: utf-8 -*-
"""
Tests of ufl to C++ expression formatting rules.

These rules can format a subset of the UFL expression language directly to C++.
The subset is currently not clearly defined, but it's basically only scalar operations
and the set of operators in UFL that coincide with C operators and cmath functions.
"""

import ufl
from ufl.constantvalue import as_ufl
from ufl.algorithms import preprocess_expression, expand_indices

from ffc.uflacs.codeutils.cpp_expr_formatting_rules import CppFormattingRules
from ffc.uflacs.codeutils.expr_formatter import ExprFormatter


def expr2cpp(expr, variables=None):
    "This is a test specific function for formatting ufl to C++."

    # Preprocessing expression before applying formatting.
    # In a compiler, one should probably assume that these
    # have been applied and use ExprFormatter directly.
    expr_data = preprocess_expression(expr)

    # Applying expand_indices instead of going through
    # compiler algorithm to get scalar valued expressions
    # that can be directlly formatted into C++
    expr = expand_indices(expr_data.preprocessed_expr)

    # This formatter is a multifunction implementing target
    # specific C++ formatting rules
    cpp_formatter = ToyCppLanguageFormatter(dh, ir)

    # This final formatter implements a generic framework handling indices etc etc.
    expr_formatter = ExprFormatter(cpp_formatter, variables or {})

    # Finally perform the formatting
    code = expr_formatter.visit(expr)
    return code


def test_cpp_formatting_of_literals():
    # Test literals
    assert expr2cpp(ufl.as_ufl(2)) == "2"
    assert expr2cpp(ufl.as_ufl(3.14)) == '3.14'
    assert expr2cpp(ufl.as_ufl(0)) == "0"
    # These are actually converted to int before formatting:
    assert expr2cpp(ufl.Identity(2)[0, 0]) == "1"
    assert expr2cpp(ufl.Identity(2)[0, 1]) == "0"
    assert expr2cpp(ufl.Identity(2)[1, 0]) == "0"
    assert expr2cpp(ufl.Identity(2)[1, 1]) == "1"
    assert expr2cpp(ufl.PermutationSymbol(3)[1, 2, 3]) == "1"
    assert expr2cpp(ufl.PermutationSymbol(3)[2, 1, 3]) == "-1"
    assert expr2cpp(ufl.PermutationSymbol(3)[1, 1, 3]) == "0"


def test_cpp_formatting_of_geometry():
    # Test geometry quantities (faked for testing!)
    x = ufl.SpatialCoordinate(ufl.interval)[0]
    assert expr2cpp(x) == "x[0]"
    x, y = ufl.SpatialCoordinate(ufl.triangle)
    assert expr2cpp(x) == "x[0]"
    assert expr2cpp(y) == "x[1]"
    nx, ny = ufl.FacetNormal(ufl.triangle)
    assert expr2cpp(nx) == "n[0]"
    assert expr2cpp(ny) == "n[1]"
    Kv = ufl.CellVolume(ufl.triangle)
    assert expr2cpp(Kv) == "volume"
    Kh = ufl.CellDiameter(ufl.triangle)
    assert expr2cpp(Kh) == "diameter"
    Kr = ufl.Circumradius(ufl.triangle)
    assert expr2cpp(Kr) == "circumradius"


def test_cpp_formatting_of_form_arguments():
    # Test form arguments (faked for testing!)
    V = ufl.FiniteElement("CG", ufl.triangle, 1)
    f = ufl.Coefficient(V, count=0)
    assert expr2cpp(f) == "w0"
    v = ufl.TestFunction(V)
    assert expr2cpp(v) == "v0"

    V = ufl.VectorElement("CG", ufl.triangle, 1)
    f = ufl.Coefficient(V, count=1)
    assert expr2cpp(f[0]) == "w0[0]"  # Renumbered to 0...
    v = ufl.Argument(V, number=3)
    assert expr2cpp(v[1]) == "v3[1]"  # NOT renumbered to 0...

    V = ufl.TensorElement("CG", ufl.triangle, 1)
    f = ufl.Coefficient(V, count=2)
    assert expr2cpp(f[1, 0]) == "w0[1][0]"  # Renumbered to 0...
    v = ufl.Argument(V, number=3)
    assert expr2cpp(v[0, 1]) == "v3[0][1]"  # NOT renumbered to 0...

    # TODO: Test mixed functions
    # TODO: Test tensor functions with symmetries


def test_cpp_formatting_of_arithmetic():
    x, y = ufl.SpatialCoordinate(ufl.triangle)
    # Test basic arithmetic operators
    assert expr2cpp(x + 3) == "3 + x[0]"
    assert expr2cpp(x * 2) == "2 * x[0]"
    assert expr2cpp(x / 2) == "x[0] / 2"
    assert expr2cpp(x * x) == "pow(x[0], 2)"  # TODO: Will gcc optimize this to x*x for us?
    assert expr2cpp(x**3) == "pow(x[0], 3)"
    # TODO: Test all basic operators


def test_cpp_formatting_of_cmath():
    x, y = ufl.SpatialCoordinate(ufl.triangle)
    # Test cmath functions
    assert expr2cpp(ufl.exp(x)) == "exp(x[0])"
    assert expr2cpp(ufl.ln(x)) == "log(x[0])"
    assert expr2cpp(ufl.sqrt(x)) == "sqrt(x[0])"
    assert expr2cpp(abs(x)) == "fabs(x[0])"
    assert expr2cpp(ufl.sin(x)) == "sin(x[0])"
    assert expr2cpp(ufl.cos(x)) == "cos(x[0])"
    assert expr2cpp(ufl.tan(x)) == "tan(x[0])"
    assert expr2cpp(ufl.asin(x)) == "asin(x[0])"
    assert expr2cpp(ufl.acos(x)) == "acos(x[0])"
    assert expr2cpp(ufl.atan(x)) == "atan(x[0])"


def test_cpp_formatting_of_derivatives():
    xx = ufl.SpatialCoordinate(ufl.triangle)
    x, y = xx
    # Test derivatives of basic operators
    assert expr2cpp(x.dx(0)) == "1"
    assert expr2cpp(x.dx(1)) == "0"
    assert expr2cpp(ufl.grad(xx)[0, 0]) == "1"
    assert expr2cpp(ufl.grad(xx)[0, 1]) == "0"
    assert expr2cpp(ufl.sin(x).dx(0)) == "cos(x[0])"

    # Test derivatives of target specific test fakes
    V = ufl.FiniteElement("CG", ufl.triangle, 1)
    f = ufl.Coefficient(V, count=0)
    assert expr2cpp(f.dx(0)) == "d1_w0[0]"
    v = ufl.Argument(V, number=3)
    assert expr2cpp(v.dx(1)) == "d1_v3[1]"  # NOT renumbered to 0...
    # TODO: Test more derivatives
    # TODO: Test variable derivatives using diff


def test_cpp_formatting_of_conditionals():
    x, y = ufl.SpatialCoordinate(ufl.triangle)
    # Test conditional expressions
    assert expr2cpp(ufl.conditional(ufl.lt(x, 2), y, 3)) \
        == "x[0] < 2 ? x[1]: 3"
    assert expr2cpp(ufl.conditional(ufl.gt(x, 2), 4 + y, 3)) \
        == "x[0] > 2 ? 4 + x[1]: 3"
    assert expr2cpp(ufl.conditional(ufl.And(ufl.le(x, 2), ufl.ge(y, 4)), 7, 8)) \
        == "x[0] <= 2 && x[1] >= 4 ? 7: 8"
    assert expr2cpp(ufl.conditional(ufl.Or(ufl.eq(x, 2), ufl.ne(y, 4)), 7, 8)) \
        == "x[0] == 2 || x[1] != 4 ? 7: 8"
    # TODO: Some tests of nested conditionals with correct precedences?


def test_cpp_formatting_precedence_handling():
    x, y = ufl.SpatialCoordinate(ufl.triangle)
    # Test precedence handling with sums
    # Note that the automatic sorting is reflected in formatting!
    assert expr2cpp(y + (2 + x)) == "x[1] + (2 + x[0])"
    assert expr2cpp((x + 2) + y) == "x[1] + (2 + x[0])"

    assert expr2cpp((2 + x) + (3 + y)) == "(2 + x[0]) + (3 + x[1])"

    assert expr2cpp((x + 3) + 2 + y) == "x[1] + (2 + (3 + x[0]))"
    assert expr2cpp(2 + (x + 3) + y) == "x[1] + (2 + (3 + x[0]))"
    assert expr2cpp(2 + (3 + x) + y) == "x[1] + (2 + (3 + x[0]))"
    assert expr2cpp(y + (2 + (3 + x))) == "x[1] + (2 + (3 + x[0]))"

    assert expr2cpp(2 + x + 3 + y) == "x[1] + (3 + (2 + x[0]))"
    assert expr2cpp(2 + x + 3 + y) == "x[1] + (3 + (2 + x[0]))"

    # Test precedence handling with divisions
    # This is more stable than sums since there is no sorting.
    assert expr2cpp((x / 2) / 3) == "(x[0] / 2) / 3"
    assert expr2cpp(x / (y / 3)) == "x[0] / (x[1] / 3)"
    assert expr2cpp((x / 2) / (y / 3)) == "(x[0] / 2) / (x[1] / 3)"
    assert expr2cpp(x / (2 / y) / 3) == "(x[0] / (2 / x[1])) / 3"

    # Test precedence handling with highest level types
    assert expr2cpp(ufl.sin(x)) == "sin(x[0])"
    assert expr2cpp(ufl.cos(x + 2)) == "cos(2 + x[0])"
    assert expr2cpp(ufl.tan(x / 2)) == "tan(x[0] / 2)"
    assert expr2cpp(ufl.acos(x + 3 * y)) == "acos(x[0] + 3 * x[1])"
    assert expr2cpp(ufl.asin(ufl.atan(x**4))) == "asin(atan(pow(x[0], 4)))"
    assert expr2cpp(ufl.sin(y) + ufl.tan(x)) == "sin(x[1]) + tan(x[0])"

    # Test precedence handling with mixed types
    assert expr2cpp(3 * (2 + x)) == "3 * (2 + x[0])"
    assert expr2cpp((2 * x) + (3 * y)) == "2 * x[0] + 3 * x[1]"
    assert expr2cpp(2 * (x + 3) * y) == "x[1] * (2 * (3 + x[0]))"
    assert expr2cpp(2 * (x + 3)**4 * y) == "x[1] * (2 * pow(3 + x[0], 4))"
    # TODO: More tests covering all types and more combinations!


def test_cpp_formatting_with_variables():
    x, y = ufl.SpatialCoordinate(ufl.triangle)
    # Test user-provided C variables for subexpressions
    # we can use variables for x[0], and sum, and power
    assert expr2cpp(x**2 + y**2, variables={x**2: 'x2', y**2: 'y2'}) == "x2 + y2"
    assert expr2cpp(x**2 + y**2, variables={x: 'z', y**2: 'y2'}) == "pow(z, 2) + y2"
    assert expr2cpp(x**2 + y**2, variables={x**2 + y**2: 'q'}) == "q"
    # we can use variables in conditionals
    assert expr2cpp(ufl.conditional(ufl.Or(ufl.eq(x, 2), ufl.ne(y, 4)), 7, 8),
                    variables={ufl.eq(x, 2): 'c1', ufl.ne(y, 4): 'c2'}) == "c1 || c2 ? 7: 8"
    # we can replace coefficients (formatted by user provided code)
    V = ufl.FiniteElement("CG", ufl.triangle, 1)
    f = ufl.Coefficient(V, count=0)
    assert expr2cpp(f, variables={f: 'f'}) == "f"
    assert expr2cpp(f**3, variables={f: 'f'}) == "pow(f, 3)"
    # variables do not replace derivatives of variable expressions
    assert expr2cpp(f.dx(0), variables={f: 'f'}) == "d1_w0[0]"

    # This test depends on which differentiation algorithm is in use
    # in UFL, representing derivatives as SpatialDerivative or Grad:
    # variables do replace variable expressions that are themselves derivatives
    assert expr2cpp(f.dx(0), variables={f.dx(0): 'df', ufl.grad(f)[0]: 'df'}) == "df"
    assert expr2cpp(ufl.grad(f)[0], variables={f.dx(0): 'df', ufl.grad(f)[0]: 'df'}) == "df"

    # TODO: Test variables in more situations with indices and derivatives

# TODO: Test various compound operators
