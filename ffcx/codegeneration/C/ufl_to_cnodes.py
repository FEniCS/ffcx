# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Tools for C/C++ expression formatting."""

import ffcx.codegeneration.lnodes as Language
import logging

import ufl

logger = logging.getLogger("ffcx")

# Table of handled math functions for different scalar types

math_table = {
    "double": {
        "sqrt": "sqrt",
        "abs": "fabs",
        "cos": "cos",
        "sin": "sin",
        "tan": "tan",
        "acos": "acos",
        "asin": "asin",
        "atan": "atan",
        "cosh": "cosh",
        "sinh": "sinh",
        "tanh": "tanh",
        "acosh": "acosh",
        "asinh": "asinh",
        "atanh": "atanh",
        "power": "pow",
        "exp": "exp",
        "ln": "log",
        "erf": "erf",
        "atan_2": "atan2",
        "min_value": "fmin",
        "max_value": "fmax",
    },
    "float": {
        "sqrt": "sqrtf",
        "abs": "fabsf",
        "cos": "cosf",
        "sin": "sinf",
        "tan": "tanf",
        "acos": "acosf",
        "asin": "asinf",
        "atan": "atanf",
        "cosh": "coshf",
        "sinh": "sinhf",
        "tanh": "tanhf",
        "acosh": "acoshf",
        "asinh": "asinhf",
        "atanh": "atanhf",
        "power": "powf",
        "exp": "expf",
        "ln": "logf",
        "erf": "erff",
        "atan_2": "atan2f",
        "min_value": "fminf",
        "max_value": "fmaxf",
    },
    "long double": {
        "sqrt": "sqrtl",
        "abs": "fabsl",
        "cos": "cosl",
        "sin": "sinl",
        "tan": "tanl",
        "acos": "acosl",
        "asin": "asinl",
        "atan": "atanl",
        "cosh": "coshl",
        "sinh": "sinhl",
        "tanh": "tanhl",
        "acosh": "acoshl",
        "asinh": "asinhl",
        "atanh": "atanhl",
        "power": "powl",
        "exp": "expl",
        "ln": "logl",
        "erf": "erfl",
        "atan_2": "atan2l",
        "min_value": "fminl",
        "max_value": "fmaxl",
    },
    "double _Complex": {
        "sqrt": "csqrt",
        "abs": "cabs",
        "cos": "ccos",
        "sin": "csin",
        "tan": "ctan",
        "acos": "cacos",
        "asin": "casin",
        "atan": "catan",
        "cosh": "ccosh",
        "sinh": "csinh",
        "tanh": "ctanh",
        "acosh": "cacosh",
        "asinh": "casinh",
        "atanh": "catanh",
        "power": "cpow",
        "exp": "cexp",
        "ln": "clog",
        "real": "creal",
        "imag": "cimag",
        "conj": "conj",
        "max_value": "fmax",
        "min_value": "fmin",
    },
    "float _Complex": {
        "sqrt": "csqrtf",
        "abs": "cabsf",
        "cos": "ccosf",
        "sin": "csinf",
        "tan": "ctanf",
        "acos": "cacosf",
        "asin": "casinf",
        "atan": "catanf",
        "cosh": "ccoshf",
        "sinh": "csinhf",
        "tanh": "ctanhf",
        "acosh": "cacoshf",
        "asinh": "casinhf",
        "atanh": "catanhf",
        "power": "cpowf",
        "exp": "cexpf",
        "ln": "clogf",
        "real": "crealf",
        "imag": "cimagf",
        "conj": "conjf",
        "max_value": "fmaxf",
        "min_value": "fminf",
    },
}


class UFL2CNodesTranslatorCpp(object):
    """UFL to CNodes translator class."""

    def __init__(self, scalar_type="double"):
        self.force_floats = False
        self.enable_strength_reduction = False
        self.scalar_type = scalar_type

        # Lookup table for handler to call when the "get" method (below) is
        # called, depending on the first argument type.
        self.call_lookup = {
            ufl.constantvalue.IntValue: self.int_value,
            ufl.constantvalue.FloatValue: self.float_value,
            ufl.constantvalue.ComplexValue: self.complex_value,
            ufl.constantvalue.Zero: self.zero,
            ufl.algebra.Product: self.product,
            ufl.algebra.Sum: self.sum,
            ufl.algebra.Division: self.division,
            ufl.algebra.Abs: self.math_function,
            ufl.algebra.Power: self.math_function,
            ufl.algebra.Real: self.math_function,
            ufl.algebra.Imag: self.math_function,
            ufl.algebra.Conj: self.math_function,
            ufl.classes.GT: self.gt,
            ufl.classes.GE: self.ge,
            ufl.classes.EQ: self.eq,
            ufl.classes.NE: self.ne,
            ufl.classes.LT: self.lt,
            ufl.classes.LE: self.le,
            ufl.classes.AndCondition: self.and_condition,
            ufl.classes.OrCondition: self.or_condition,
            ufl.classes.NotCondition: self.not_condition,
            ufl.classes.Conditional: self.conditional,
            ufl.classes.MinValue: self.math_function,
            ufl.classes.MaxValue: self.math_function,
            ufl.mathfunctions.Sqrt: self.math_function,
            ufl.mathfunctions.Ln: self.math_function,
            ufl.mathfunctions.Exp: self.math_function,
            ufl.mathfunctions.Cos: self.math_function,
            ufl.mathfunctions.Sin: self.math_function,
            ufl.mathfunctions.Tan: self.math_function,
            ufl.mathfunctions.Cosh: self.math_function,
            ufl.mathfunctions.Sinh: self.math_function,
            ufl.mathfunctions.Tanh: self.math_function,
            ufl.mathfunctions.Acos: self.math_function,
            ufl.mathfunctions.Asin: self.math_function,
            ufl.mathfunctions.Atan: self.math_function,
            ufl.mathfunctions.Erf: self.math_function,
            ufl.mathfunctions.Atan2: self.math_function,
            ufl.mathfunctions.MathFunction: self.math_function,
            ufl.mathfunctions.BesselJ: self.bessel_j,
            ufl.mathfunctions.BesselY: self.bessel_y,
        }

    def get(self, o, *args):
        # Call appropriate handler, depending on the type of o
        otype = type(o)
        if otype in self.call_lookup:
            return self.call_lookup[otype](o, *args)
        else:
            raise RuntimeError(f"Missing C formatting rule for expr type {otype}.")

    def expr(self, o, *args):
        """Raise generic fallback with error message for missing rules."""
        raise RuntimeError(f"Missing C formatting rule for expr type {o._ufl_class_}.")

    # === Formatting rules for scalar literals ===

    def zero(self, o):
        return Language.LiteralFloat(0.0)

    def float_value(self, o):
        return Language.LiteralFloat(float(o))

    def int_value(self, o):
        if self.force_floats:
            return self.float_value(o)
        return Language.LiteralInt(int(o))

    def complex_value(self, o):
        return Language.LiteralFloat(o.value())

    # === Formatting rules for arithmetic operators ===

    def sum(self, o, a, b):
        return Language.Add(a, b)

    def product(self, o, a, b):
        return Language.Product([a, b])

    def division(self, o, a, b):
        if self.enable_strength_reduction:
            return Language.Product([a, Language.Div(1.0, b)])
        else:
            return Language.Div(a, b)

    # === Formatting rules for conditional expressions ===

    def conditional(self, o, c, t, f):
        return Language.Conditional(c, t, f)

    def eq(self, o, a, b):
        return Language.EQ(a, b)

    def ne(self, o, a, b):
        return Language.NE(a, b)

    def le(self, o, a, b):
        return Language.LE(a, b)

    def ge(self, o, a, b):
        return Language.GE(a, b)

    def lt(self, o, a, b):
        return Language.LT(a, b)

    def gt(self, o, a, b):
        return Language.GT(a, b)

    def and_condition(self, o, a, b):
        return Language.And(a, b)

    def or_condition(self, o, a, b):
        return Language.Or(a, b)

    def not_condition(self, o, a):
        return Language.Not(a)

    # === Formatting rules for cmath functions ===

    def math_function(self, o, *args):
        # Fallback for unhandled MathFunction subclass:
        # attempting to just call it.
        return Language.MathFunction(o._ufl_handler_name_, args)

    # def math_function(self, o, *args):
    #     k = o._ufl_handler_name_
    #     try:
    #         name = math_table[self.scalar_type].get(k)
    #     except Exception as e:
    #         raise type(e)("Math function not found:", self.scalar_type, k)
    #     if name is None:
    #         raise RuntimeError("Not supported in current scalar mode")
    #     return Language.MathFunction(name, args)

    # === Formatting rules for bessel functions ===
    # Some Bessel functions exist in gcc, as XSI extensions
    # but not all.

    def bessel_j(self, o, n, v):
        assert "complex" not in self.scalar_type
        n = int(float(n))
        if n == 0:
            return Language.Call("j0", v)
        elif n == 1:
            return Language.Call("j1", v)
        else:
            return Language.Call("jn", (n, v))

    def bessel_y(self, o, n, v):
        assert "complex" not in self.scalar_type
        n = int(float(n))
        if n == 0:
            return Language.Call("y0", v)
        elif n == 1:
            return Language.Call("y1", v)
        else:
            return Language.Call("yn", (n, v))
