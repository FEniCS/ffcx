# -*- coding: utf-8 -*-
# Copyright (C) 2011-2017 Martin Sandve Alnæs
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

from ufl.corealg.multifunction import MultiFunction
#from ufl.corealg.map_dag import map_expr_dag


class UFL2CNodesMixin(object):
    """Rules collection mixin for a UFL to CNodes translator class."""
    def __init__(self, language):
        self.L = language

        self.force_floats = False
        self.enable_strength_reduction = False

    # === Error handlers for missing formatting rules ===

    def expr(self, o):
        "Generic fallback with error message for missing rules."
        error("Missing C++ formatting rule for expr type {0}.".format(o._ufl_class_))

    # === Formatting rules for scalar literals ===

    def zero(self, o):
        return self.L.LiteralFloat(0.0)

    #def complex_value(self, o):
    #    return self.L.ComplexValue(complex(o))

    def float_value(self, o):
        return self.L.LiteralFloat(float(o))

    def int_value(self, o):
        if self.force_floats:
            return self.float_value(o)
        return self.L.LiteralInt(int(o))

    # === Formatting rules for arithmetic operators ===

    def sum(self, o, a, b):
        return self.L.Add(a, b)

    #def sub(self, o, a, b): # Not in UFL
    #    return self.L.Sub(a, b)

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
        # Fallback for unhandled MathFunction subclass: attempting to just call it.
        # TODO: Introduce a UserFunction to UFL to keep it separate from MathFunction?
        return self.L.Call(o._name, op)

    def sqrt(self, o, op):
        return self._cmath("sqrt", op)

    #def cbrt(self, o, op):  # Not in UFL
    #    return self._cmath("cbrt", op)

    # cmath also has log10 etc
    def ln(self, o, op):
        return self._cmath("log", op)

    # cmath also has exp2 etc
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

    #def acosh(self, o, op):  # Not in UFL
    #    return self._cmath("acosh", op)

    #def asinh(self, o, op):  # Not in UFL
    #    return self._cmath("asinh", op)

    #def atanh(self, o, op):  # Not in UFL
    #    return self._cmath("atanh", op)

    def erf(self, o, op):
        return self._cmath("erf", op)

    #def erfc(self, o, op):  # Not in UFL
    #    # C++11 stl has this function
    #    return self._cmath("erfc", op)


class RulesForC(object):
    def _cmath(self, name, op):
        return self.L.Call(name, op)

    def power(self, o, a, b):
        return self.L.Call("pow", (a, b))

    def abs(self, o, op):
        return self.L.Call("fabs", op)

    def min_value(self, o, a, b):
        return self.L.Call("fmin", (a, b))

    def max_value(self, o, a, b):
        return self.L.Call("fmax", (a, b))

    # ignoring bessel functions


class RulesForCpp(object):
    def _cmath(self, name, op):
        return self.L.Call("std::" + name, op)

    def power(self, o, a, b):
        return self.L.Call("std::pow", (a, b))

    def abs(self, o, op):
        return self.L.Call("std::abs", op)

    def min_value(self, o, a, b):
        return self.L.Call("std::min", (a, b))

    def max_value(self, o, a, b):
        return self.L.Call("std::max", (a, b))

    # === Formatting rules for bessel functions ===

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


class UFL2CNodesTranslatorC(MultiFunction, UFL2CNodesMixin, RulesForC):
    """UFL to CNodes translator class."""
    def __init__(self, language):
        MultiFunction.__init__(self)
        UFL2CNodesMixin.__init__(self, language)


class UFL2CNodesTranslatorCpp(MultiFunction, UFL2CNodesMixin, RulesForCpp):
    """UFL to CNodes translator class."""
    def __init__(self, language):
        MultiFunction.__init__(self)
        UFL2CNodesMixin.__init__(self, language)
