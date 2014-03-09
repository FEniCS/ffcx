
from uflacs.utils.log import uflacs_assert, error

import ufl

class CppFormattingRules(object):
    """C++ formatting rules collection.

    This is the base class for target specific cpp formatter class.
    """
    def __init__(self):
        pass

    # === Utility functions ===

    # TODO: Other solution to this
    def add_include(self, include, system):
        "Callback to be implemented by subclass for recording includes needed."
        pass

    # TODO: Other solution to this
    def add_using(self, using):
        "Callback to be implemented by subclass for recording using symbols needed."
        pass

    # === Error handlers for missing formatting rules ===

    # Generic fallback error messages for missing rules:
    def expr(self, o):
        error("Missing C++ formatting rule for expr type {0}.".format(o._uflclass))

    def terminal(self, o, mt):
        error("Missing C++ formatting rule for terminal type {0}.".format(o._uflclass))

    # Unexcepted type checks:
    def variable(self, o, *ops):
        error("Expecting variables to be removed before formatting C++ code.")
        return ops[0] # or just fall through like this if necessary

    def invalid_request(self, o, *ops):
        error("Invalid request for C++ formatting of a {0}, str = {1}".format(o._uflclass, str(o)))
    wrapper_type = invalid_request
    index_sum    = invalid_request
    indexed      = invalid_request
    derivative   = invalid_request
    restricted   = invalid_request

    argument            = invalid_request
    coefficient         = invalid_request
    geometric_quantitiy = invalid_request

    # === Formatting rules for literal constants ===

    def constant_value(self, e, mt):
        error("Missing C++ rule for constant value type {0}.".format(e._uflclass))

    def zero(self, e, mt):
        return "0"

    def int_value(self, e, mt):
        return "{0}".format(int(e))

    def float_value(self, e, mt):
        # Using configurable precision parameter from ufl
        if mt.global_derivatives or mt.local_derivatives:
            return self.zero(None)
        else:
            return ufl.constantvalue.format_float(float(e))

    def identity(self, o, mt):
        if mt.global_derivatives or mt.local_derivatives:
            return "0"
        else:
            return "1" if mt.component[0] == mt.component[1] else "0"

    # === Formatting rules for arithmetic operators ===

    def sum(self, o, *ops):
        return " + ".join(ops)

    def product(self, o, *ops):
        return " * ".join(ops)

    def division(self, o, a, b):
        return "{0} / {1}".format(a, b)

    # === Formatting rules for cmath functions ===

    def _cmath(self, name, op):
        self.add_using("std::{0}".format(name))
        return "{0}({1})".format(name, op)

    def math_function(self, o, op):
        return self._cmath(o._name, op)

    def power(self, o, a, b):
        self.add_using("std::pow")
        return "pow({0}, {1})".format(a, b)

    def sqrt(self, o, op):
        return self._cmath("sqrt", op)

    def ln(self, o, op):
        return self._cmath("log", op)

    def exp(self, o, op):
        return self._cmath("exp", op)

    def abs(self, o, op):
        return "fabs({0})".format(op)

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
        return "erf({0})".format(op)

    # === Formatting rules for bessel functions ===

    def _bessel(self, o, n, v, name):
        self.add_include("boost/math/special_functions.hpp", True)
        self.add_using("boost::math::{0}".format(name))
        return "{0}({1}, {2})".format(name, n, v)

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
