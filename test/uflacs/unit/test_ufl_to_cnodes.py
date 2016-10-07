# -*- coding: utf-8 -*-

import ufl
from ufl import *
from ufl import as_ufl

import ffc.uflacs.language
from ffc.uflacs.language.ufl_to_cnodes import UFL2CNodesTranslatorCpp, UFL2CNodesTranslatorC


def test_ufl_to_cnodes():

    L = ffc.uflacs.language.cnodes
    translate = UFL2CNodesTranslatorCpp(L)

    f = ufl.CellVolume(ufl.triangle)
    g = ufl.CellVolume(ufl.triangle)
    h = ufl.CellVolume(ufl.triangle)

    x = L.Symbol("x")
    y = L.Symbol("y")
    z = L.Symbol("z")

    examples = [
        (f + g, (x, y), "x + y"),
        # (f-g, (x,y), "x - y"), # - is not currently a UFL operator
        (f * g, (x, y), "x * y"),
        (f / g, (x, y), "x / y"),
        (f**g, (x, y), "std::pow(x, y)"),
        (exp(f), (x,), "std::exp(x)"),
        (ln(f), (x,), "std::log(x)"),
        (abs(f), (x,), "std::abs(x)"),
        (sqrt(f), (x,), "std::sqrt(x)"),
        #(cbrt(f), (x,), "std::cbrt(x)"),
        (sin(f), (x,), "std::sin(x)"),
        (cos(f), (x,), "std::cos(x)"),
        (tan(f), (x,), "std::tan(x)"),
        (asin(f), (x,), "std::asin(x)"),
        (acos(f), (x,), "std::acos(x)"),
        (atan(f), (x,), "std::atan(x)"),
        (sinh(f), (x,), "std::sinh(x)"),
        (cosh(f), (x,), "std::cosh(x)"),
        (tanh(f), (x,), "std::tanh(x)"),
        (erf(f), (x,), "std::erf(x)"),
        #(erfc(f), (x,), "std::erfc(x)"),
        (min_value(f, g), (x, y), "std::min(x, y)"),
        (max_value(f, g), (x, y), "std::max(x, y)"),
        (bessel_I(1, g), (x, y), "boost::math::cyl_bessel_i(x, y)"),
        (bessel_J(1, g), (x, y), "boost::math::cyl_bessel_j(x, y)"),
        (bessel_K(1, g), (x, y), "boost::math::cyl_bessel_k(x, y)"),
        (bessel_Y(1, g), (x, y), "boost::math::cyl_neumann(x, y)"),
        (f < g, (x, y), "x < y"),
        (f > g, (x, y), "x > y"),
        (f <= g, (x, y), "x <= y"),
        (f >= g, (x, y), "x >= y"),
        (eq(f, g), (x, y), "x == y"),
        (ne(f, g), (x, y), "x != y"),
        (And(f < g, f > g), (x, y), "x && y"),
        (Or(f < g, f > g), (x, y), "x || y"),
        (Not(f < g), (x,), "!x"),
        (conditional(f < g, g, h), (x, y, z), "x ? y : z"),
    ]
    for expr, args, code in examples:
        # Warning: This test is subtle: translate will look at the type of expr and
        #  ignore its operands, i.e. not translating the full tree but only one level.
        assert str(translate(expr, *args)) == code


    # C specific translation:
    translate = UFL2CNodesTranslatorC(L)
    examples = [
        (sin(f), (x,), "sin(x)"),
        (f**g, (x, y), "pow(x, y)"),
        (exp(f), (x,), "exp(x)"),
        (abs(f), (x,), "fabs(x)"),
        (min_value(f, g), (x, y), "fmin(x, y)"),
        (max_value(f, g), (x, y), "fmax(x, y)"),
        ]
    for expr, args, code in examples:
        # Warning: This test is subtle: translate will look at the type of expr and
        #  ignore its operands, i.e. not translating the full tree but only one level.
        assert str(translate(expr, *args)) == code
