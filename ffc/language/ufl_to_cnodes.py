# -*- coding: utf-8 -*-
# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Tools for C/C++ expression formatting."""

import logging

from ffc import FFCError
from ufl.corealg.multifunction import MultiFunction

logger = logging.getLogger(__name__)

# Table of handled math functions in real and complex modes
math_table = {'sqrt': ('sqrt', 'csqrt'),
              'abs': ('fabs', 'cabs'),
              'cos': ('cos', 'ccos'),
              'sin': ('sin', 'csin'),
              'tan': ('tan', 'ctan'),
              'acos': ('acos', 'cacos'),
              'asin': ('asin', 'casin'),
              'atan': ('atan', 'catan'),
              'cosh': ('cosh', 'ccosh'),
              'sinh': ('sinh', 'csinh'),
              'tanh': ('tanh', 'ctanh'),
              'acosh': ('acosh', 'cacosh'),
              'asinh': ('asinh', 'casinh'),
              'atanh': ('atanh', 'catanh'),
              'power': ('pow', 'cpow'),
              'exp': ('exp', 'cexp'),
              'ln': ('log', 'clog'),
              'real': (None, 'creal'),
              'imag': (None, 'cimag'),
              'conj': (None, 'conj'),
              'erf': ('erf', None),
              'atan_2': ('atan2', None),
              'min_value': ('fmin', None),
              'max_value': ('fmax', None)}


class UFL2CNodesTranslatorCpp(MultiFunction):
    """UFL to CNodes translator class."""

    def __init__(self, language, complex_mode=False):
        MultiFunction.__init__(self)

        self.L = language
        self.force_floats = False
        self.enable_strength_reduction = False
        self.complex_mode = 1 if complex_mode else 0

    # === Error handlers for missing formatting rules ===

    def expr(self, o):
        """Generic fallback with error message for missing rules."""
        raise FFCError("Missing C++ formatting rule for expr type {0}.".format(o._ufl_class_))

    # === Formatting rules for scalar literals ===

    def zero(self, o):
        return self.L.LiteralFloat(0.0)

    def float_value(self, o):
        return self.L.LiteralFloat(float(o))

    def int_value(self, o):
        if self.force_floats:
            return self.float_value(o)
        return self.L.LiteralInt(int(o))

    def complex_value(self, o):
        return self.L.LiteralFloat(o.value())

    # === Formatting rules for arithmetic operators ===

    def sum(self, o, a, b):
        return self.L.Add(a, b)

    def product(self, o, a, b):
        return self.L.Mul(a, b)

    def division(self, o, a, b):
        if self.enable_strength_reduction:
            return self.L.Mul(a, self.L.Div(1.0, b))
        else:
            return self.L.Div(a, b)

    # === Formatting rules for conditional expressions ===

    def conditional(self, o, c, t, f):
        return self.L.Conditional(c, t, f)

    def eq(self, o, a, b):
        return self.L.EQ(a, b)

    def ne(self, o, a, b):
        return self.L.NE(a, b)

    def le(self, o, a, b):
        return self.L.LE(a, b)

    def ge(self, o, a, b):
        return self.L.GE(a, b)

    def lt(self, o, a, b):
        return self.L.LT(a, b)

    def gt(self, o, a, b):
        return self.L.GT(a, b)

    def and_condition(self, o, a, b):
        return self.L.And(a, b)

    def or_condition(self, o, a, b):
        return self.L.Or(a, b)

    def not_condition(self, o, a):
        return self.L.Not(a)

    # === Formatting rules for cmath functions ===

    def math_function(self, o, op):
        # Fallback for unhandled MathFunction subclass:
        # attempting to just call it.
        return self.L.Call(o._name, op)

    def _cmath(self, o, *args):
        k = o._ufl_handler_name_
        try:
            name = math_table[k]
        except Exception as e:
            raise type(e)("Math function not found:", k)
        name = name[self.complex_mode]
        if name is None:
            raise RuntimeError("Not supported in current complex mode")
        return self.L.Call(name, args)

    real = _cmath
    imag = _cmath
    conj = _cmath
    sqrt = _cmath
    ln = _cmath
    exp = _cmath
    cos = _cmath
    sin = _cmath
    tan = _cmath
    cosh = _cmath
    sinh = _cmath
    tanh = _cmath
    acos = _cmath
    asin = _cmath
    atan = _cmath
    erf = _cmath
    power = _cmath
    abs = _cmath
    atan_2 = _cmath
    min_value = _cmath
    max_value = _cmath

    # === Formatting rules for bessel functions ===
    # Some Bessel functions exist in gcc, as XSI extensions
    # but not all.

    def bessel_j(self, o, n, v):
        assert self.complex_mode == 0
        n = int(float(n))
        if n == 0:
            return self.L.Call("j0", v)
        elif n == 1:
            return self.L.Call("j1", v)
        else:
            return self.L.Call("jn", (n, v))

    def bessel_y(self, o, n, v):
        assert self.complex_mode == 0
        n = int(float(n))
        if n == 0:
            return self.L.Call("y0", v)
        elif n == 1:
            return self.L.Call("y1", v)
        else:
            return self.L.Call("yn", (n, v))
