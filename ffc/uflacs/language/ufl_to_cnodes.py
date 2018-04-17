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


class UFL2CNodesTranslatorCpp(MultiFunction):
    """UFL to CNodes translator class."""

    def __init__(self, language):
        MultiFunction.__init__(self)

        self.L = language
        self.force_floats = False
        self.enable_strength_reduction = False

    # === Error handlers for missing formatting rules ===

    def expr(self, o):
        "Generic fallback with error message for missing rules."
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

    def _cmath(self, name, op):
        return self.L.Call(name, op)

    def sqrt(self, o, op):
        return self._cmath("sqrt", op)

    def ln(self, o, op):
        return self._cmath("log", op)

    def exp(self, o, op):
        return self._cmath("exp", op)

    def cos(self, o, op):
        return self._cmath("cos", op)

    def sin(self, o, op):
        return self._cmath("sin", op)

    def tan(self, o, op):
        return self._cmath("tan", op)

    def cosh(self, o, op):
        return self._cmath("cosh", op)

    def sinh(self, o, op):
        return self._cmath("sinh", op)

    def tanh(self, o, op):
        return self._cmath("tanh", op)

    def atan_2(self, o, y, x):
        return self._cmath("atan2", (y, x))

    def acos(self, o, op):
        return self._cmath("acos", op)

    def asin(self, o, op):
        return self._cmath("asin", op)

    def atan(self, o, op):
        return self._cmath("atan", op)

    def erf(self, o, op):
        return self._cmath("erf", op)

    def power(self, o, a, b):
        return self._cmath("pow", (a, b))

    def abs(self, o, op):
        return self._cmath("fabs", op)

    def min_value(self, o, a, b):
        return self._cmath("min", (a, b))

    def max_value(self, o, a, b):
        return self._cmath("max", (a, b))

    # === Formatting rules for bessel functions ===
    # Currently in boost; will be in C++17

    def _bessel(self, o, n, v, name):
        return self.L.Call("boost::math::" + name, (n, v))

    def bessel_i(self, o, n, v):
        return self._bessel(o, n, v, "cyl_bessel_i")

    def bessel_j(self, o, n, v):
        return self._bessel(o, n, v, "cyl_bessel_j")

    def bessel_k(self, o, n, v):
        return self._bessel(o, n, v, "cyl_bessel_k")

    def bessel_y(self, o, n, v):
        return self._bessel(o, n, v, "cyl_neumann")
