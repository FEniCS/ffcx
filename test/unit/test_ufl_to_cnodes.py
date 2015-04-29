
import ufl
from ufl import *
from ufl import as_ufl

import uflacs.language
from uflacs.language.ufl_to_cnodes import UFL2CNodesTranslator

def test_ufl_to_cnodes():

    L = uflacs.language.cnodes
    translate = UFL2CNodesTranslator(L)

    f = ufl.CellVolume(ufl.triangle)
    g = ufl.CellVolume(ufl.triangle)
    h = ufl.CellVolume(ufl.triangle)

    x = L.Symbol("x")
    y = L.Symbol("y")
    z = L.Symbol("z")

    examples = [
        (f+g, (x,y), "x + y"),
        #(f-g, (x,y), "x - y"), # - is not currently a UFL operator
        (f*g, (x,y), "x * y"),
        (f/g, (x,y), "x / y"),
        (f**g, (x,y), "pow(x, y)"),
        (exp(f), (x,), "exp(x)"),
        (ln(f), (x,), "log(x)"),
        (abs(f), (x,), "abs(x)"),
        (sqrt(f), (x,), "sqrt(x)"),
        #(cbrt(f), (x,), "cbrt(x)"),
        (sin(f), (x,), "sin(x)"),
        (cos(f), (x,), "cos(x)"),
        (tan(f), (x,), "tan(x)"),
        (asin(f), (x,), "asin(x)"),
        (acos(f), (x,), "acos(x)"),
        (atan(f), (x,), "atan(x)"),
        (sinh(f), (x,), "sinh(x)"),
        (cosh(f), (x,), "cosh(x)"),
        (tanh(f), (x,), "tanh(x)"),
        (erf(f), (x,), "erf(x)"),
        #(erfc(f), (x,), "erfc(x)"),
        (min_value(f,g), (x,y), "min(x, y)"),
        (max_value(f,g), (x,y), "max(x, y)"),
        (bessel_I(1, g), (x,y), "cyl_bessel_i(x, y)"),
        (bessel_J(1, g), (x,y), "cyl_bessel_j(x, y)"),
        (bessel_K(1, g), (x,y), "cyl_bessel_k(x, y)"),
        (bessel_Y(1, g), (x,y), "cyl_neumann(x, y)"),
        (f < g, (x,y), "x < y"),
        (f > g, (x,y), "x > y"),
        (f <= g, (x,y), "x <= y"),
        (f >= g, (x,y), "x >= y"),
        (eq(f, g), (x,y), "x == y"),
        (ne(f, g), (x,y), "x != y"),
        (And(f<g, f>g), (x,y), "x && y"),
        (Or(f<g, f>g), (x,y), "x || y"),
        (Not(f<g), (x,), "!x"),
        (conditional(f<g,g,h), (x,y,z), "x ? y : z"),
        ]
    for expr, args, code in examples:
        # Warning: This test is subtle: translate will look at the type of expr and
        #  ignore its operands, i.e. not translating the full tree but only one level.
        assert str(translate(expr, *args)) == code
