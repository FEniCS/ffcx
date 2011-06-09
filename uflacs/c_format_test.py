
from log import info, warning, error
from output import uflacs_assert

import ufl
from ufl.algorithms.transformations import MultiFunction, Transformer
from ufl.algorithms import expand_derivatives, expand_indices

class CppTestFormatter(object):
    """Example cpp formatter class, used for the test cases.
    Override the same functions for your particular target."""

    def spatial_coordinate(self, o, component=(), derivatives=()):
        if len(derivatives) > 1:
            return "0"
        i = component[0] if component else 0
        if derivatives:
            d, = derivatives
            return "1" if i == d else "0"
        else:
            return "x[%d]" % i

    def facet_normal(self, o, component=(), derivatives=()):
        if derivatives:
            return "0"
        i = component[0] if component else 0
        return "n[%d]" % i

    def cell_volume(self, o, component=(), derivatives=()):
        uflacs_assert(not component, "Expecting no component for scalar value.")
        if derivatives:
            return "0"
        return "K_vol"

    def circumradius(self, o, component=(), derivatives=()):
        uflacs_assert(not component, "Expecting no component for scalar value.")
        if derivatives:
            return "0"
        return "K_rad"

    def coefficient(self, o, component=(), derivatives=()):
        basename = "dw" if derivatives else "w"
        code = "%s[%d]" % (basename, o.count())
        if component:
            code += " " + "".join("[%d]" % c for c in component)
        if derivatives:
            uflacs_assert(len(derivatives) == 1, "Not implemented.")
            code += " " + "".join("[%d]" % d for d in derivatives)
        return code

    def argument(self, o, component=(), derivatives=()):
        basename = "dv" if derivatives else "v"
        code = "%s[%d]" % (basename, o.count())
        if component:
            code += " " + "".join("[%d]" % c for c in component)
        if derivatives:
            uflacs_assert(len(derivatives) == 1, "Not implemented.")
            code += " " + "".join("[%d]" % d for d in derivatives)
        return code

def format_expression_as_cpp(expr, target_formatter=None, variables=None):
    "This is a test specific function for formatting ufl to C++."
    from c_format import CppFormatterRules, CodeFormatter

    # This formatter is a multifunction with single operator
    # formatting rules for generic C++ formatting
    cpp_formatter = CppFormatterRules(target_formatter)

    # This final formatter implements a generic framework handling indices etc etc.
    variables = variables or {}
    code_formatter = CodeFormatter(cpp_formatter, variables)
    code = code_formatter.visit(expr)
    return code

def format_expression_as_test_cpp(expr, variables=None):
    "This is a test specific function for formatting ufl to C++."

    # Preprocessing expression before applying formatting.
    # In a compiler, one should probably assume that these
    # have been applied and use CodeFormatter directly.
    e = expr
    e = expand_derivatives(e)
    e = expand_indices(e)

    # This formatter is a multifunction implementing target
    # specific formatting rules, here mocked for testing.
    test_formatter = CppTestFormatter()

    # Call the generic function
    return format_expression_as_cpp(e, test_formatter, variables)

class TestStats:
    def __init__(self):
        self.fails = 0
        self.successes = 0

    def __str__(self):
        return "Fails: %d, Successes: %d." % (self.fails, self.successes)

def test_code_formatting():
    teststats = TestStats()
    def assertEqual(expr, code, variables=None):
        r = format_expression_as_test_cpp(expr, variables)
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

    # Test geometry quantities (faked for testing!)
    x = ufl.cell1D.x
    assertEqual(x, "x[0]")
    x, y = ufl.cell2D.x
    assertEqual(x, "x[0]")
    assertEqual(y, "x[1]")
    nx, ny = ufl.cell2D.n
    assertEqual(nx, "n[0]")
    assertEqual(ny, "n[1]")
    Kv = ufl.cell2D.volume
    assertEqual(Kv, "K_vol")
    Kr = ufl.cell2D.circumradius
    assertEqual(Kr, "K_rad")

    # Test form arguments (faked for testing!)
    V = ufl.FiniteElement("CG", ufl.cell2D, 1)
    f = ufl.Coefficient(V).reconstruct(count=0)
    assertEqual(f, "w[0]")
    v = ufl.Argument(V).reconstruct(count=0)
    assertEqual(v, "v[0]")

    V = ufl.VectorElement("CG", ufl.cell2D, 1)
    f = ufl.Coefficient(V).reconstruct(count=1)
    assertEqual(f[0], "w[1] [0]")
    v = ufl.Argument(V).reconstruct(count=3)
    assertEqual(v[1], "v[3] [1]")

    V = ufl.TensorElement("CG", ufl.cell2D, 1)
    f = ufl.Coefficient(V).reconstruct(count=2)
    assertEqual(f[1,0], "w[2] [1][0]")
    v = ufl.Argument(V).reconstruct(count=3)
    assertEqual(v[0,1], "v[3] [0][1]")

    # TODO: Test mixed functions
    # TODO: Test tensor functions with symmetries

    # Test basic arithmetic operators
    assertEqual(x + 3, "3 + x[0]")
    assertEqual(x * 2, "2 * x[0]")
    assertEqual(x / 2, "x[0] / 2")
    assertEqual(x*x, "std::pow(x[0], 2)") # TODO: Will gcc optimize this to x*x for us?
    assertEqual(x**3, "std::pow(x[0], 3)")
    # TODO: Test all basic operators

    # Test cmath functions
    assertEqual(ufl.exp(x), "std::exp(x[0])")
    assertEqual(ufl.ln(x), "std::log(x[0])")
    assertEqual(ufl.sqrt(x), "std::sqrt(x[0])")
    assertEqual(abs(x), "std::abs(x[0])")
    assertEqual(ufl.sin(x), "std::sin(x[0])")
    assertEqual(ufl.cos(x), "std::cos(x[0])")
    assertEqual(ufl.tan(x), "std::tan(x[0])")
    assertEqual(ufl.asin(x), "std::asin(x[0])")
    assertEqual(ufl.acos(x), "std::acos(x[0])")
    assertEqual(ufl.atan(x), "std::atan(x[0])")

    # Test derivatives of basic operators
    assertEqual(ufl.Identity(2)[0,0].dx(0), "0")
    assertEqual(x.dx(0), "1")
    assertEqual(x.dx(1), "0")
    assertEqual(ufl.sin(x).dx(0), "std::cos(x[0])")

    # Test derivatives of target specific test fakes
    V = ufl.FiniteElement("CG", ufl.cell2D, 1)
    f = ufl.Coefficient(V).reconstruct(count=0)
    assertEqual(f.dx(0), "dw[0] [0]")
    v = ufl.Argument(V).reconstruct(count=0)
    assertEqual(v.dx(0), "dv[0] [0]")
    # TODO: Test more derivatives
    # TODO: Test variable derivatives using diff

    # Test conditional expressions
    assertEqual(ufl.conditional(ufl.lt(x, 2), y, 3),
                "x[0] < 2 ? x[1]: 3")
    assertEqual(ufl.conditional(ufl.gt(x, 2), 4+y, 3),
                "x[0] > 2 ? 4 + x[1]: 3")
    assertEqual(ufl.conditional(ufl.And(ufl.le(x, 2), ufl.ge(y, 4)), 7, 8),
                "x[0] <= 2 && x[1] >= 4 ? 7: 8")
    assertEqual(ufl.conditional(ufl.Or(ufl.eq(x, 2), ufl.ne(y, 4)), 7, 8),
                "x[0] == 2 || x[1] != 4 ? 7: 8")
    # TODO: Some tests of nested conditionals with correct precedences?

    # Test precedence handling with sums
    # Note that the automatic sorting is reflected in formatting!
    assertEqual(y + (2 + x), "x[1] + (2 + x[0])")
    assertEqual((x + 2) + y, "x[1] + (2 + x[0])")

    assertEqual((2 + x) + (3 + y), "(2 + x[0]) + (3 + x[1])")

    assertEqual((x + 3) + 2 + y, "x[1] + (2 + (3 + x[0]))")
    assertEqual(2 + (x + 3) + y, "x[1] + (2 + (3 + x[0]))")
    assertEqual(2 + (3 + x) + y, "x[1] + (2 + (3 + x[0]))")
    assertEqual(y + (2 + (3 + x)), "x[1] + (2 + (3 + x[0]))")

    assertEqual(2 + x + 3 + y, "x[1] + (3 + (2 + x[0]))")
    assertEqual(2 + x + 3 + y, "x[1] + (3 + (2 + x[0]))")

    # Test precedence handling with divisions
    # This is more stable than sums since there is no sorting.
    assertEqual((x / 2) / 3, "(x[0] / 2) / 3")
    assertEqual(x / (y / 3), "x[0] / (x[1] / 3)")
    assertEqual((x / 2) / (y / 3), "(x[0] / 2) / (x[1] / 3)")
    assertEqual(x / (2 / y) / 3, "(x[0] / (2 / x[1])) / 3")

    # Test precedence handling with highest level types
    assertEqual(ufl.sin(x), "std::sin(x[0])")
    assertEqual(ufl.cos(x+2), "std::cos(2 + x[0])")
    assertEqual(ufl.tan(x/2), "std::tan(x[0] / 2)")
    assertEqual(ufl.acos(x + 3 * y), "std::acos(x[0] + 3 * x[1])")
    assertEqual(ufl.asin(ufl.atan(x**4)), "std::asin(std::atan(std::pow(x[0], 4)))")
    assertEqual(ufl.sin(y) + ufl.tan(x), "std::sin(x[1]) + std::tan(x[0])")

    # Test precedence handling with mixed types
    assertEqual(3 * (2 + x), "3 * (2 + x[0])")
    assertEqual((2 * x) + (3 * y), "2 * x[0] + 3 * x[1]")
    assertEqual(2 * (x + 3) * y, "x[1] * (2 * (3 + x[0]))")
    assertEqual(2 * (x + 3)**4 * y, "x[1] * (2 * std::pow(3 + x[0], 4))")
    # TODO: More tests covering all types and more combinations!

    # Test user-provided C variables for subexpressions
    # we can use variables for x[0], and sum, and power
    assertEqual(x**2 + y**2, "x2 + y2", variables={x**2: 'x2', y**2: 'y2'})
    assertEqual(x**2 + y**2, "std::pow(z, 2) + y2", variables={x: 'z', y**2: 'y2'})
    assertEqual(x**2 + y**2, "q", variables={x**2 + y**2: 'q'})
    # we can use variables in conditionals
    assertEqual(ufl.conditional(ufl.Or(ufl.eq(x, 2), ufl.ne(y, 4)), 7, 8),
                "c1 || c2 ? 7: 8",
                variables={ufl.eq(x, 2): 'c1', ufl.ne(y, 4): 'c2'})
    # we can replace coefficients (formatted by user provided code)
    V = ufl.FiniteElement("CG", ufl.cell2D, 1)
    f = ufl.Coefficient(V).reconstruct(count=0)
    assertEqual(f, "f", variables={f: 'f'})
    assertEqual(f**3, "std::pow(f, 3)", variables={f: 'f'})
    # variables do not replace derivatives of variable expressions
    assertEqual(f.dx(0), "dw[0] [0]", variables={f: 'f'})
    # variables do replace variable expressions that are themselves derivatives
    assertEqual(f.dx(0), "df", variables={f.dx(0): 'df'})

    # TODO: Test variables in more situations with indices and derivatives

    # TODO: Test various compound operators

    # Report final test results
    print teststats
    return teststats.fails
