# -*- coding: utf-8 -*-
"""
Tests of LaTeX formatting rules.
"""

from ffc.log import ffc_assert, info, warning, error
from ffc.uflacs.codeutils.expr_formatter import ExprFormatter
from ffc.uflacs.codeutils.latex_expr_formatting_rules import LatexFormatter

import ufl
from ufl.algorithms import preprocess_expression, expand_indices


def expr2latex(expr, variables=None):
    "This is a test specific function for formatting ufl to LaTeX."

    # Preprocessing expression before applying formatting.  In a
    # compiler, one should probably assume that these have been
    # applied and use ExprFormatter directly.
    expr_data = preprocess_expression(expr)
    expr = expand_indices(expr_data.preprocessed_expr)

    # This formatter is a multifunction with single operator
    # formatting rules for generic LaTeX formatting
    latex_formatter = LatexFormatter()

    # This final formatter implements a generic framework handling
    # indices etc etc.
    variables = variables or {}
    expr_formatter = ExprFormatter(latex_formatter, variables)
    code = expr_formatter.visit(expr)
    return code


def assertLatexEqual(self, expr, code, variables=None):
    r = expr2latex(expr, variables)
    self.assertEqual(code, r)


def test_latex_formatting_of_literals():
    # Test literals
    assert expr2latex(ufl.as_ufl(2)) == "2"
    assert expr2latex(ufl.as_ufl(3.14)) == '3.14'
    assert expr2latex(ufl.as_ufl(0)) == "0"
    # These are actually converted to int before formatting:
    assert expr2latex(ufl.Identity(2)[0, 0]) == "1"
    assert expr2latex(ufl.Identity(2)[0, 1]) == "0"
    assert expr2latex(ufl.Identity(2)[1, 0]) == "0"
    assert expr2latex(ufl.Identity(2)[1, 1]) == "1"
    assert expr2latex(ufl.PermutationSymbol(3)[1, 2, 3]) == "1"
    assert expr2latex(ufl.PermutationSymbol(3)[2, 1, 3]) == "-1"
    assert expr2latex(ufl.PermutationSymbol(3)[1, 1, 3]) == "0"


def test_latex_formatting_of_geometry():
    # Test geometry quantities
    x = ufl.SpatialCoordinate(ufl.interval)[0]
    assert expr2latex(x) == "x_0"
    x, y = ufl.SpatialCoordinate(ufl.triangle)
    assert expr2latex(x) == "x_0"
    assert expr2latex(y) == "x_1"
    nx, ny = ufl.FacetNormal(ufl.triangle)
    assert expr2latex(nx) == "n_0"
    assert expr2latex(ny) == "n_1"
    Kv = ufl.CellVolume(ufl.triangle)
    assert expr2latex(Kv) == r"K_{\text{vol}}"
    Kr = ufl.Circumradius(ufl.triangle)
    assert expr2latex(Kr) == r"K_{\text{rad}}"


def test_latex_formatting_of_form_arguments():
    # Test form arguments (faked for testing!)
    U = ufl.FiniteElement("CG", ufl.triangle, 1)
    V = ufl.VectorElement("CG", ufl.triangle, 1)
    W = ufl.TensorElement("CG", ufl.triangle, 1)

    v = ufl.TestFunction(U)
    assert expr2latex(v) == r"\overset{0}{v}"
    f = ufl.Coefficient(U, count=0)
    assert expr2latex(f) == r"\overset{0}{w}"

    f = ufl.Coefficient(V, count=1)
    assert expr2latex(f[0]) == r"\overset{1}{w}_{0}"  # NOT renumbered to 0...
    v = ufl.Argument(V, number=3)
    assert expr2latex(v[1]) == r"\overset{3}{v}_{1}"  # NOT renumbered to 0...

    f = ufl.Coefficient(W, count=2)
    assert expr2latex(f[1, 0]) == r"\overset{2}{w}_{1 0}"  # NOT renumbered to 0...
    v = ufl.Argument(W, number=3)
    assert expr2latex(v[0, 1]) == r"\overset{3}{v}_{0 1}"  # NOT renumbered to 0...

    # TODO: Test mixed functions
    # TODO: Test tensor functions with symmetries


def test_latex_formatting_of_arithmetic():
    x = ufl.SpatialCoordinate(ufl.triangle)[0]
    assert expr2latex(x + 3) == "3 + x_0"
    assert expr2latex(x * 2) == "2 x_0"
    assert expr2latex(x / 2) == r"\frac{x_0}{2}"
    assert expr2latex(x * x) == r"{x_0}^{2}"  # TODO: Will gcc optimize this to x*x for us?
    assert expr2latex(x**3) == r"{x_0}^{3}"
    # TODO: Test all basic operators


def test_latex_formatting_of_cmath():
    x = ufl.SpatialCoordinate(ufl.triangle)[0]
    assert expr2latex(ufl.exp(x)) == r"e^{x_0}"
    assert expr2latex(ufl.ln(x)) == r"\ln(x_0)"
    assert expr2latex(ufl.sqrt(x)) == r"\sqrt{x_0}"
    assert expr2latex(abs(x)) == r"\|x_0\|"
    assert expr2latex(ufl.sin(x)) == r"\sin(x_0)"
    assert expr2latex(ufl.cos(x)) == r"\cos(x_0)"
    assert expr2latex(ufl.tan(x)) == r"\tan(x_0)"
    assert expr2latex(ufl.asin(x)) == r"\arcsin(x_0)"
    assert expr2latex(ufl.acos(x)) == r"\arccos(x_0)"
    assert expr2latex(ufl.atan(x)) == r"\arctan(x_0)"


def test_latex_formatting_of_derivatives():
    xx = ufl.SpatialCoordinate(ufl.triangle)
    x = xx[0]
    # Test derivatives of basic operators
    assert expr2latex(x.dx(0)) == "1"
    assert expr2latex(x.dx(1)) == "0"
    assert expr2latex(ufl.grad(xx)[0, 0]) == "1"
    assert expr2latex(ufl.grad(xx)[0, 1]) == "0"
    assert expr2latex(ufl.sin(x).dx(0)) == r"\cos(x_0)"

    # Test derivatives of form arguments
    V = ufl.FiniteElement("CG", ufl.triangle, 1)
    f = ufl.Coefficient(V, count=0)
    assert expr2latex(f.dx(0)) == r"\overset{0}{w}_{, 0}"
    v = ufl.Argument(V, number=3)
    assert expr2latex(v.dx(1)) == r"\overset{3}{v}_{, 1}"
    # TODO: Test more derivatives
    # TODO: Test variable derivatives using diff


def xtest_latex_formatting_of_conditionals():
    # Test conditional expressions
    assert expr2latex(ufl.conditional(ufl.lt(x, 2), y, 3)) == "x_0 < 2 ? x_1: 3"
    assert expr2latex(ufl.conditional(ufl.gt(x, 2), 4 + y, 3)) == "x_0 > 2 ? 4 + x_1: 3"
    assert expr2latex(ufl.conditional(ufl.And(ufl.le(x, 2), ufl.ge(y, 4)), 7, 8)) == "x_0 <= 2 && x_1 >= 4 ? 7: 8"
    assert expr2latex(ufl.conditional(ufl.Or(ufl.eq(x, 2), ufl.ne(y, 4)), 7, 8)) == "x_0 == 2 || x_1 != 4 ? 7: 8"
    # TODO: Some tests of nested conditionals with correct precedences?


def test_latex_formatting_precedence_handling():
    x, y = ufl.SpatialCoordinate(ufl.triangle)
    # Test precedence handling with sums
    # Note that the automatic sorting is reflected in formatting!
    assert expr2latex(y + (2 + x)) == "x_1 + (2 + x_0)"
    assert expr2latex((x + 2) + y) == "x_1 + (2 + x_0)"

    assert expr2latex((2 + x) + (3 + y)) == "(2 + x_0) + (3 + x_1)"

    assert expr2latex((x + 3) + 2 + y) == "x_1 + (2 + (3 + x_0))"
    assert expr2latex(2 + (x + 3) + y) == "x_1 + (2 + (3 + x_0))"
    assert expr2latex(2 + (3 + x) + y) == "x_1 + (2 + (3 + x_0))"
    assert expr2latex(y + (2 + (3 + x))) == "x_1 + (2 + (3 + x_0))"

    assert expr2latex(2 + x + 3 + y) == "x_1 + (3 + (2 + x_0))"
    assert expr2latex(2 + x + 3 + y) == "x_1 + (3 + (2 + x_0))"

    # Test precedence handling with divisions
    # This is more stable than sums since there is no sorting.
    assert expr2latex((x / 2) / 3) == r"\frac{(\frac{x_0}{2})}{3}"
    assert expr2latex(x / (y / 3)) == r"\frac{x_0}{(\frac{x_1}{3})}"
    assert expr2latex((x / 2) / (y / 3)) == r"\frac{(\frac{x_0}{2})}{(\frac{x_1}{3})}"
    assert expr2latex(x / (2 / y) / 3) == r"\frac{(\frac{x_0}{(\frac{2}{x_1})})}{3}"

    # Test precedence handling with highest level types
    assert expr2latex(ufl.sin(x)) == r"\sin(x_0)"
    assert expr2latex(ufl.cos(x + 2)) == r"\cos(2 + x_0)"
    assert expr2latex(ufl.tan(x / 2)) == r"\tan(\frac{x_0}{2})"
    assert expr2latex(ufl.acos(x + 3 * y)) == r"\arccos(x_0 + 3 x_1)"
    assert expr2latex(ufl.asin(ufl.atan(x**4))) == r"\arcsin(\arctan({x_0}^{4}))"
    assert expr2latex(ufl.sin(y) + ufl.tan(x)) == r"\sin(x_1) + \tan(x_0)"

    # Test precedence handling with mixed types
    assert expr2latex(3 * (2 + x)) == "3 (2 + x_0)"
    assert expr2latex((2 * x) + (3 * y)) == "2 x_0 + 3 x_1"
    assert expr2latex(2 * (x + 3) * y) == "x_1 (2 (3 + x_0))"


def _fixme():
    assert expr2latex(2 * (x + 3)**4 * y) == "x_1 (2 {(3 + x_0)}^{4)"  # FIXME: Precedence handling fails here
    # TODO: More tests covering all types and more combinations!


def test_latex_formatting_of_variables():
    x, y = ufl.SpatialCoordinate(ufl.triangle)
    # Test user-provided C variables for subexpressions
    # we can use variables for x[0], and sum, and power
    assert expr2latex(x**2 + y**2, variables={x**2: 'x2', y**2: 'y2'}) == "x2 + y2"
    assert expr2latex(x**2 + y**2, variables={x: 'z', y**2: 'y2'}) == r"{z}^{2} + y2"
    assert expr2latex(x**2 + y**2, variables={x**2 + y**2: 'q'}) == "q"
    # we can use variables in conditionals
    if 0:
        assert expr2latex(ufl.conditional(ufl.Or(ufl.eq(x, 2), ufl.ne(y, 4)), 7, 8), variables={ufl.eq(x, 2): 'c1', ufl.ne(y, 4): 'c2'}) == "c1 || c2 ? 7: 8"
    # we can replace coefficients (formatted by user provided code)
    V = ufl.FiniteElement("CG", ufl.triangle, 1)
    f = ufl.Coefficient(V, count=0)
    assert expr2latex(f, variables={f: 'f'}) == "f"
    assert expr2latex(f**3, variables={f: 'f'}) == r"{f}^{3}"
    # variables do not replace derivatives of variable expressions
    assert expr2latex(f.dx(0), variables={f: 'f'}) == r"\overset{0}{w}_{, 0}"
    # variables do replace variable expressions that are themselves derivatives
    assert expr2latex(f.dx(0), variables={f.dx(0): 'df', ufl.grad(f)[0]: 'df'}) == "df"
    assert expr2latex(ufl.grad(f)[0], variables={f.dx(0): 'df', ufl.grad(f)[0]: 'df'}) == "df"

    # TODO: Test variables in more situations with indices and derivatives

# TODO: Test various compound operators
