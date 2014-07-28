
from ffc.log import error

import ufl


class CppFormattingRules(object):

    """C++ formatting rules collection.

    This is the base class for target specific cpp formatter class.
    """

    def __init__(self):
        pass

    # === Error handlers for missing formatting rules ===

    # Generic fallback error messages for missing rules:
    def expr(self, o):
        error("Missing C++ formatting rule for expr type {0}.".format(o._ufl_class_))

    def terminal(self, o, mt):
        error("Missing C++ formatting rule for terminal type {0}.".format(o._ufl_class_))

    # Unexcepted type checks:
    def variable(self, o, *ops):
        error("Expecting variables to be removed before formatting C++ code.")
        return ops[0]  # or just fall through like this if necessary

    def invalid_request(self, o, *ops):
        error("Invalid request for C++ formatting of a {0}, str = {1}".format(o._ufl_class_, str(o)))
    wrapper_type = invalid_request
    index_sum = invalid_request
    indexed = invalid_request
    derivative = invalid_request
    restricted = invalid_request

    argument = invalid_request
    coefficient = invalid_request
    geometric_quantitiy = invalid_request

    # === Formatting rules for literal constants ===

    def constant_value(self, e, mt=None):
        error("Missing C++ rule for constant value type {0}.".format(e._ufl_class_))

    def zero(self, e, mt=None):
        return "0.0"

    def int_value(self, e, mt=None):
        if mt is not None and (mt.global_derivatives or mt.local_derivatives):
            return "0"
        else:
            return "{0}".format(int(e))

    def float_value(self, e, mt=None):
        if mt is not None and (mt.global_derivatives or mt.local_derivatives):
            return "0.0"
        else:
            # Using configurable precision parameter from ufl
            return ufl.constantvalue.format_float(float(e))

    def identity(self, o, mt):
        if mt.global_derivatives or mt.local_derivatives:
            return "0.0"
        else:
            return "1.0" if mt.component[0] == mt.component[1] else "0.0"

    # === Formatting rules for arithmetic operators ===

    def sum(self, o, *ops):
        return " + ".join(ops)

    def product(self, o, *ops):
        return " * ".join(ops)

    def division(self, o, a, b):
        return "{0} / {1}".format(a, b)

    # === Formatting rules for cmath functions ===

    def _cmath(self, name, op):
        return "{0}{1}({2})".format("std::", name, op)

    def math_function(self, o, op):
        return self._cmath(o._name, op)

    def power(self, o, a, b):
        return "{0}pow({1}, {2})".format("std::", a, b)

    def sqrt(self, o, op):
        return self._cmath("sqrt", op)

    def ln(self, o, op):
        return self._cmath("log", op)

    def exp(self, o, op):
        return self._cmath("exp", op)

    def abs(self, o, op):
        # return "fabs({0})".format(op) # C version
        return self._cmath("abs", op)  # C++ stl version

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
        # return self._cmath("erf", op) # C++11 stl has this function
        return "erf({0})".format(op)

    # def erfc(self, o, op): # Not in UFL
    # return self._cmath("erfc", op) # C++11 stl has this function
    #    return "erfc({0})".format(op)

    def min_value(self, o, a, b):
        # return "fmin({0}, {1})".format(a, b) # C99 version
        return "{0}min({1}, {2})".format("std::", a, b)  # C++ stl version

    def max_value(self, o, a, b):
        # return "fmax({0}, {1})".format(a, b) # C99 version
        return "{0}max({1}, {2})".format("std::", a, b)  # C++ stl version

    # === Formatting rules for bessel functions ===

    def _bessel(self, o, n, v, name):
        return "{0}{1}({2}, {3})".format("boost::math::", name, n, v)

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
        return "{0} ? {1}: {2}".format(c, t, f)

    def eq(self, o, a, b):
        return "{0} == {1}".format(a, b)

    def ne(self, o, a, b):
        return "{0} != {1}".format(a, b)

    def le(self, o, a, b):
        return "{0} <= {1}".format(a, b)

    def ge(self, o, a, b):
        return "{0} >= {1}".format(a, b)

    def lt(self, o, a, b):
        return "{0} < {1}".format(a, b)

    def gt(self, o, a, b):
        return "{0} > {1}".format(a, b)

    def and_condition(self, o, a, b):
        return "{0} && {1}".format(a, b)

    def or_condition(self, o, a, b):
        return "{0} || {1}".format(a, b)

    def not_condition(self, o, a):
        return "!{0}".format(a)

from ufl.algorithms import MultiFunction


class CppExprFormatter(MultiFunction, CppFormattingRules):

    """C++ formatter class. Customize by copying this class and overriding specific rules."""

    def __init__(self):
        MultiFunction.__init__(self)
        CppFormattingRules.__init__(self)
