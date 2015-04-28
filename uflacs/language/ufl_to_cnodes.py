# -*- coding: utf-8 -*-
# Copyright (C) 2011-2015 Martin Sandve Alnæs
#
# This file is part of UFLACS.
#
# UFLACS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# UFLACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with UFLACS. If not, see <http://www.gnu.org/licenses/>.

"""Tools for C/C++ expression formatting."""

from ffc.log import error

import ufl

from ufl.corealg.multifunction import MultiFunction
from ufl.corealg.map_dag import map_expr_dag

class UFL2CNodesMixin(object):
    """Rules collection mixin for a UFL to CNodes translator class."""
    def __init__(self, language):
        self.lang = language

    # === Error handlers for missing formatting rules ===

    # Generic fallback error messages for missing rules:
    def expr(self, o):
        error("Missing C++ formatting rule for expr type {0}.".format(o._ufl_class_))

    # === Formatting rules for arithmetic operators ===

    def sum(self, o, a, b):
        return self.lang.Add(a, b)

    def product(self, o, a, b):
        return self.lang.Mul(a, b)

    def division(self, o, a, b):
        return self.lang.Div(a, b)

    # === Formatting rules for cmath functions ===

    def power(self, o, a, b):
        return self.lang.Call("pow", (a, b))

    def _cmath(self, name, op):
        # TODO: Configurable namespacing?
        #name = "std::" + name
        return self.lang.Call(name, op)

    def math_function(self, o, op):
        return self._cmath(o._name, op)

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

    def acos(self, o, op):
        return self._cmath("acos", op)

    def asin(self, o, op):
        return self._cmath("asin", op)

    def atan(self, o, op):
        return self._cmath("atan", op)

    def erf(self, o, op):
        #return self._cmath("erf", op) # C++11 stl has this function
        return self.lang.Call("erf", op)

    #def erfc(self, o, op): # Not in UFL
    #    return self._cmath("erfc", op) # C++11 stl has this function

    # FIXME: Need to incorporate some functions as CNode types,
    #        to switch between C and C++ behaviour:
    def abs(self, o, op):
        #return Call("fabs", op) # C version
        return self._cmath("abs", op)  # C++ stl version

    def min_value(self, o, a, b):
        #return self.lang.Call("fmin", (a, b)) # C99 version
        return self.lang.Call("min", (a, b))  # C++ stl version

    def max_value(self, o, a, b):
        #return self.lang.Call("fmax", (a, b)) # C99 version
        return self.lang.Call("max", (a, b))  # C++ stl version

    # === Formatting rules for bessel functions ===

    def _bessel(self, o, n, v, name):
        #return "{0}{1}({2}, {3})".format("boost::math::", name, n, v)
        #name = "boost::math::" + name
        return self.lang.Call(name, (n, v))

    def bessel_i(self, o, n, v):
        return self._bessel(o, n, v, "cyl_bessel_i")

    def bessel_j(self, o, n, v):
        return self._bessel(o, n, v, "cyl_bessel_j")

    def bessel_k(self, o, n, v):
        return self._bessel(o, n, v, "cyl_bessel_k")

    def bessel_y(self, o, n, v):
        return self._bessel(o, n, v, "cyl_neumann")

    # === Formatting rules for conditional expressions ===

    def conditional(self, o, c, t, f):
        return self.lang.Conditional(c, t, f)

    def eq(self, o, a, b):
        return self.lang.EQ(a, b)

    def ne(self, o, a, b):
        return self.lang.NE(a, b)

    def le(self, o, a, b):
        return self.lang.LE(a, b)

    def ge(self, o, a, b):
        return self.lang.GE(a, b)

    def lt(self, o, a, b):
        return self.lang.LT(a, b)

    def gt(self, o, a, b):
        return self.lang.GT(a, b)

    def and_condition(self, o, a, b):
        return self.lang.And(a, b)

    def or_condition(self, o, a, b):
        return self.lang.Or(a, b)

    def not_condition(self, o, a):
        return self.lang.Not(a)


class UFL2CNodesTranslator(MultiFunction, UFL2CNodesMixin):
    """UFL to CNodes translator class.

    Customize by copying this class and overriding specific rules."""
    def __init__(self, language):
        MultiFunction.__init__(self)
        UFL2CNodesMixin.__init__(self, language)


def ufl_to_cnodes(expr):
    import uflacs.language.cnodes # FIXME: Move to uflacs
    language = uflacs.language.cnodes # FIXME: Make argument
    return map_expr_dag(UFL2CNodesTranslator(language), expr)
