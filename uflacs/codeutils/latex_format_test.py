
from uflacs.utils.log import info, warning, error
from uflacs.utils.assertions import uflacs_assert

import ufl
from ufl.algorithms import preprocess_expression

def format_expression_as_test_latex(expr, variables=None):
    "This is a test specific function for formatting ufl to LaTeX."
    from uflacs.codeutils.code_formatter import CodeFormatter
    from uflacs.codeutils.latex_format import LatexFormatterRules

    # Preprocessing expression before applying formatting.
    # In a compiler, one should probably assume that these
    # have been applied and use CodeFormatter directly.
    expr_data = preprocess_expression(expr)

    # This formatter is a multifunction with single operator
    # formatting rules for generic LaTeX formatting
    latex_formatter = LatexFormatterRules()

    # This final formatter implements a generic framework handling indices etc etc.
    variables = variables or {}
    code_formatter = CodeFormatter(latex_formatter, variables)
    code = code_formatter.visit(expr_data.preprocessed_expr)
    return code

class TestStats:
    def __init__(self):
        self.fails = 0
        self.successes = 0

    def __str__(self):
        return "Fails: %d, Successes: %d." % (self.fails, self.successes)

def test_latex_formatting():
    teststats = TestStats()
    def assertEqual(expr, code, variables=None):
        r = format_expression_as_test_latex(expr, variables)
        if code == r:
            teststats.successes += 1
        else:
            teststats.fails += 1
            print "FAIL assertEqual:"
            print "Actual:   '%s'" % r
            print "Expected: '%s'" % code
            print

    # Test literals
    assertEqual(ufl.as_ufl(2), "2")
    assertEqual(ufl.as_ufl(3.14), '3.14')
    assertEqual(ufl.as_ufl(0), "0")
    # These are actually converted to int before formatting:
    assertEqual(ufl.Identity(2)[0,0], "1")
    assertEqual(ufl.Identity(2)[0,1], "0")
    assertEqual(ufl.Identity(2)[1,0], "0")
    assertEqual(ufl.Identity(2)[1,1], "1")
    assertEqual(ufl.PermutationSymbol(3)[1,2,3], "1")
    assertEqual(ufl.PermutationSymbol(3)[2,1,3], "-1")
    assertEqual(ufl.PermutationSymbol(3)[1,1,3], "0")

    # Test geometry quantities
    x = ufl.cell1D.x
    assertEqual(x, "x_0")
    x, y = ufl.cell2D.x
    assertEqual(x, "x_0")
    assertEqual(y, "x_1")
    nx, ny = ufl.cell2D.n
    assertEqual(nx, "n_0")
    assertEqual(ny, "n_1")
    Kv = ufl.cell2D.volume
    assertEqual(Kv, r"K_{\text{vol}}")
    Kr = ufl.cell2D.circumradius
    assertEqual(Kr, r"K_{\text{rad}}")

    # Test form arguments (faked for testing!)
    V = ufl.FiniteElement("CG", ufl.cell2D, 1)
    f = ufl.Coefficient(V).reconstruct(count=0)
    assertEqual(f, r"\overset{0}{w}")
    v = ufl.Argument(V).reconstruct(count=0)
    assertEqual(v, r"\overset{0}{v}")

    V = ufl.VectorElement("CG", ufl.cell2D, 1)
    f = ufl.Coefficient(V).reconstruct(count=1)
    assertEqual(f[0], r"\overset{0}{w}_{0}") # Renumbered to 0...
    v = ufl.Argument(V).reconstruct(count=3)
    assertEqual(v[1], r"\overset{0}{v}_{1}") # Renumbered to 0...

    V = ufl.TensorElement("CG", ufl.cell2D, 1)
    f = ufl.Coefficient(V).reconstruct(count=2)
    assertEqual(f[1,0], r"\overset{0}{w}_{1 0}") # Renumbered to 0...
    v = ufl.Argument(V).reconstruct(count=3)
    assertEqual(v[0,1], r"\overset{0}{v}_{0 1}") # Renumbered to 0...

    # TODO: Test mixed functions
    # TODO: Test tensor functions with symmetries

    # Test basic arithmetic operators
    assertEqual(x + 3, "3 + x_0")
    assertEqual(x * 2, "2 x_0")
    assertEqual(x / 2, r"\frac{x_0}{2}")
    assertEqual(x*x, r"pow(x_0, 2)") # TODO: Will gcc optimize this to x*x for us?
    assertEqual(x**3, r"pow(x_0, 3)")
    # TODO: Test all basic operators

    # Test cmath functions
    assertEqual(ufl.exp(x), r"exp(x_0)")
    assertEqual(ufl.ln(x), r"\ln(x_0)")
    assertEqual(ufl.sqrt(x), r"sqrt(x_0)")
    assertEqual(abs(x), r"\|x_0\|")
    assertEqual(ufl.sin(x), r"\sin(x_0)")
    assertEqual(ufl.cos(x), r"\cos(x_0)")
    assertEqual(ufl.tan(x), r"\tan(x_0)")
    assertEqual(ufl.asin(x), r"\arcsin(x_0)")
    assertEqual(ufl.acos(x), r"\arccos(x_0)")
    assertEqual(ufl.atan(x), r"\arctan(x_0)")

    # Test derivatives of basic operators
    assertEqual(ufl.Identity(2)[0,0].dx(0), "0")
    assertEqual(x.dx(0), "1")
    assertEqual(x.dx(1), "0")
    assertEqual(ufl.sin(x).dx(0), r"\cos(x_0)")

    # Test derivatives of target specific test fakes
    V = ufl.FiniteElement("CG", ufl.cell2D, 1)
    f = ufl.Coefficient(V).reconstruct(count=0)
    assertEqual(f.dx(0), r"\overset{0}{w}_{, 0}")
    v = ufl.Argument(V).reconstruct(count=3)
    assertEqual(v.dx(1), r"\overset{0}{v}_{, 1}")
    # TODO: Test more derivatives
    # TODO: Test variable derivatives using diff

    # Test conditional expressions
    if 0:
        assertEqual(ufl.conditional(ufl.lt(x, 2), y, 3),
                    "x_0 < 2 ? x_1: 3")
        assertEqual(ufl.conditional(ufl.gt(x, 2), 4+y, 3),
                    "x_0 > 2 ? 4 + x_1: 3")
        assertEqual(ufl.conditional(ufl.And(ufl.le(x, 2), ufl.ge(y, 4)), 7, 8),
                    "x_0 <= 2 && x_1 >= 4 ? 7: 8")
        assertEqual(ufl.conditional(ufl.Or(ufl.eq(x, 2), ufl.ne(y, 4)), 7, 8),
                    "x_0 == 2 || x_1 != 4 ? 7: 8")
        # TODO: Some tests of nested conditionals with correct precedences?

    # Test precedence handling with sums
    # Note that the automatic sorting is reflected in formatting!
    assertEqual(y + (2 + x), "x_1 + (2 + x_0)")
    assertEqual((x + 2) + y, "x_1 + (2 + x_0)")

    assertEqual((2 + x) + (3 + y), "(2 + x_0) + (3 + x_1)")

    assertEqual((x + 3) + 2 + y, "x_1 + (2 + (3 + x_0))")
    assertEqual(2 + (x + 3) + y, "x_1 + (2 + (3 + x_0))")
    assertEqual(2 + (3 + x) + y, "x_1 + (2 + (3 + x_0))")
    assertEqual(y + (2 + (3 + x)), "x_1 + (2 + (3 + x_0))")

    assertEqual(2 + x + 3 + y, "x_1 + (3 + (2 + x_0))")
    assertEqual(2 + x + 3 + y, "x_1 + (3 + (2 + x_0))")

    # Test precedence handling with divisions
    # This is more stable than sums since there is no sorting.
    assertEqual((x / 2) / 3, r"\frac{(\frac{x_0}{2})}{3}")
    assertEqual(x / (y / 3), r"\frac{x_0}{(\frac{x_1}{3})}")
    assertEqual((x / 2) / (y / 3), r"\frac{(\frac{x_0}{2})}{(\frac{x_1}{3})}")
    assertEqual(x / (2 / y) / 3, r"\frac{(\frac{x_0}{(\frac{2}{x_1})})}{3}")

    # Test precedence handling with highest level types
    assertEqual(ufl.sin(x), r"\sin(x_0)")
    assertEqual(ufl.cos(x+2), r"\cos(2 + x_0)")
    assertEqual(ufl.tan(x/2), r"\tan(\frac{x_0}{2})")
    assertEqual(ufl.acos(x + 3 * y), r"\arccos(x_0 + 3 x_1)")
    assertEqual(ufl.asin(ufl.atan(x**4)), r"\arcsin(\arctan(pow(x_0, 4)))")
    assertEqual(ufl.sin(y) + ufl.tan(x), r"\sin(x_1) + \tan(x_0)")

    # Test precedence handling with mixed types
    assertEqual(3 * (2 + x), "3 (2 + x_0)")
    assertEqual((2 * x) + (3 * y), "2 x_0 + 3 x_1")
    assertEqual(2 * (x + 3) * y, "x_1 (2 (3 + x_0))")
    assertEqual(2 * (x + 3)**4 * y, "x_1 (2 pow(3 + x_0, 4))")
    # TODO: More tests covering all types and more combinations!

    # Test user-provided C variables for subexpressions
    # we can use variables for x[0], and sum, and power
    assertEqual(x**2 + y**2, "x2 + y2", variables={x**2: 'x2', y**2: 'y2'})
    assertEqual(x**2 + y**2, r"pow(z, 2) + y2", variables={x: 'z', y**2: 'y2'})
    assertEqual(x**2 + y**2, "q", variables={x**2 + y**2: 'q'})
    # we can use variables in conditionals
    if 0:
        assertEqual(ufl.conditional(ufl.Or(ufl.eq(x, 2), ufl.ne(y, 4)), 7, 8),
                    "c1 || c2 ? 7: 8",
                    variables={ufl.eq(x, 2): 'c1', ufl.ne(y, 4): 'c2'})
    # we can replace coefficients (formatted by user provided code)
    V = ufl.FiniteElement("CG", ufl.cell2D, 1)
    f = ufl.Coefficient(V).reconstruct(count=0)
    assertEqual(f, "f", variables={f: 'f'})
    assertEqual(f**3, r"pow(f, 3)", variables={f: 'f'})
    # variables do not replace derivatives of variable expressions
    assertEqual(f.dx(0), r"\overset{0}{w}_{, 0}", variables={f: 'f'})
    # variables do replace variable expressions that are themselves derivatives
    assertEqual(f.dx(0), "df", variables={f.dx(0): 'df'})

    # TODO: Test variables in more situations with indices and derivatives

    # TODO: Test various compound operators

    # Report final test results
    print teststats
    return teststats.fails
