# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Tools for C/C++ expression formatting."""

import logging

import ufl

logger = logging.getLogger("ffcx")

# Table of handled math functions for different scalar types

math_table = {'double': {'sqrt': 'sqrt',
                         'abs': 'fabs',
                         'cos': 'cos',
                         'sin': 'sin',
                         'tan': 'tan',
                         'acos': 'acos',
                         'asin': 'asin',
                         'atan': 'atan',
                         'cosh': 'cosh',
                         'sinh': 'sinh',
                         'tanh': 'tanh',
                         'acosh': 'acosh',
                         'asinh': 'asinh',
                         'atanh': 'atanh',
                         'power': 'pow',
                         'exp': 'exp',
                         'ln': 'log',
                         'erf': 'erf',
                         'atan_2': 'atan2',
                         'min_value': 'fmin',
                         'max_value': 'fmax'},

              'float': {'sqrt': 'sqrtf',
                        'abs': 'fabsf',
                        'cos': 'cosf',
                        'sin': 'sinf',
                        'tan': 'tanf',
                        'acos': 'acosf',
                        'asin': 'asinf',
                        'atan': 'atanf',
                        'cosh': 'coshf',
                        'sinh': 'sinhf',
                        'tanh': 'tanhf',
                        'acosh': 'acoshf',
                        'asinh': 'asinhf',
                        'atanh': 'atanhf',
                        'power': 'powf',
                        'exp': 'expf',
                        'ln': 'logf',
                        'erf': 'erff',
                        'atan_2': 'atan2f',
                        'min_value': 'fminf',
                        'max_value': 'fmaxf'},

              'long double': {'sqrt': 'sqrtl',
                              'abs': 'fabsl',
                              'cos': 'cosl',
                              'sin': 'sinl',
                              'tan': 'tanl',
                              'acos': 'acosl',
                              'asin': 'asinl',
                              'atan': 'atanl',
                              'cosh': 'coshl',
                              'sinh': 'sinhl',
                              'tanh': 'tanhl',
                              'acosh': 'acoshl',
                              'asinh': 'asinhl',
                              'atanh': 'atanhl',
                              'power': 'powl',
                              'exp': 'expl',
                              'ln': 'logl',
                              'erf': 'erfl',
                              'atan_2': 'atan2l',
                              'min_value': 'fminl',
                              'max_value': 'fmaxl'},

              'double _Complex': {'sqrt': 'csqrt',
                                  'abs': 'cabs',
                                  'cos': 'ccos',
                                  'sin': 'csin',
                                  'tan': 'ctan',
                                  'acos': 'cacos',
                                  'asin': 'casin',
                                  'atan': 'catan',
                                  'cosh': 'ccosh',
                                  'sinh': 'csinh',
                                  'tanh': 'ctanh',
                                  'acosh': 'cacosh',
                                  'asinh': 'casinh',
                                  'atanh': 'catanh',
                                  'power': 'cpow',
                                  'exp': 'cexp',
                                  'ln': 'clog',
                                  'real': 'creal',
                                  'imag': 'cimag',
                                  'conj': 'conj',
                                  'max_value': 'fmax',
                                  'min_value': 'fmin'},

              'float _Complex': {'sqrt': 'csqrtf',
                                 'abs': 'cabsf',
                                 'cos': 'ccosf',
                                 'sin': 'csinf',
                                 'tan': 'ctanf',
                                 'acos': 'cacosf',
                                 'asin': 'casinf',
                                 'atan': 'catanf',
                                 'cosh': 'ccoshf',
                                 'sinh': 'csinhf',
                                 'tanh': 'ctanhf',
                                 'acosh': 'cacoshf',
                                 'asinh': 'casinhf',
                                 'atanh': 'catanhf',
                                 'power': 'cpowf',
                                 'exp': 'cexpf',
                                 'ln': 'clogf',
                                 'real': 'crealf',
                                 'imag': 'cimagf',
                                 'conj': 'conjf',
                                 'max_value': 'fmaxf',
                                 'min_value': 'fminf'}}


class UFL2CNodesTranslatorCpp(object):
    """UFL to CNodes translator class."""

    def __init__(self, language, scalar_type="double"):
        self.L = language
        self.force_floats = False
        self.enable_strength_reduction = False
        self.scalar_type = scalar_type

        # Lookup table for handler to call when the "get" method (below) is
        # called, depending on the first argument type.
        self.call_lookup = {ufl.constantvalue.IntValue: self.int_value,
                            ufl.constantvalue.FloatValue: self.float_value,
                            ufl.constantvalue.ComplexValue: self.complex_value,
                            ufl.constantvalue.Zero: self.zero,
                            ufl.algebra.Product: self.product,
                            ufl.algebra.Sum: self.sum,
                            ufl.algebra.Division: self.division,
                            ufl.algebra.Abs: self._cmath,
                            ufl.algebra.Power: self._cmath,
                            ufl.algebra.Real: self._cmath,
                            ufl.algebra.Imag: self._cmath,
                            ufl.algebra.Conj: self._cmath,
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
                            ufl.classes.MinValue: self._cmath,
                            ufl.classes.MaxValue: self._cmath,
                            ufl.mathfunctions.Sqrt: self._cmath,
                            ufl.mathfunctions.Ln: self._cmath,
                            ufl.mathfunctions.Exp: self._cmath,
                            ufl.mathfunctions.Cos: self._cmath,
                            ufl.mathfunctions.Sin: self._cmath,
                            ufl.mathfunctions.Tan: self._cmath,
                            ufl.mathfunctions.Cosh: self._cmath,
                            ufl.mathfunctions.Sinh: self._cmath,
                            ufl.mathfunctions.Tanh: self._cmath,
                            ufl.mathfunctions.Acos: self._cmath,
                            ufl.mathfunctions.Asin: self._cmath,
                            ufl.mathfunctions.Atan: self._cmath,
                            ufl.mathfunctions.Erf: self._cmath,
                            ufl.mathfunctions.Atan2: self._cmath,
                            ufl.mathfunctions.MathFunction: self.math_function,
                            ufl.mathfunctions.BesselJ: self.bessel_j,
                            ufl.mathfunctions.BesselY: self.bessel_y}

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
            name = math_table[self.scalar_type].get(k)
        except Exception as e:
            raise type(e)("Math function not found:", self.scalar_type, k)
        if name is None:
            raise RuntimeError("Not supported in current scalar mode")
        return self.L.Call(name, args)

    # === Formatting rules for bessel functions ===
    # Some Bessel functions exist in gcc, as XSI extensions
    # but not all.

    def bessel_j(self, o, n, v):
        assert "complex" not in self.scalar_type
        n = int(float(n))
        if n == 0:
            return self.L.Call("j0", v)
        elif n == 1:
            return self.L.Call("j1", v)
        else:
            return self.L.Call("jn", (n, v))

    def bessel_y(self, o, n, v):
        assert "complex" not in self.scalar_type
        n = int(float(n))
        if n == 0:
            return self.L.Call("y0", v)
        elif n == 1:
            return self.L.Call("y1", v)
        else:
            return self.L.Call("yn", (n, v))
