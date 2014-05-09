
from uflacs.utils.log import uflacs_assert, info, warning, error

import ufl

# TODO: Assuming in this code that preprocessed expressions
# are formatted, so no compounds etc. are included here.
# Would be nice to format e.g. dot(u, v) -> u \cdot v.

class LatexFormattingRules(object):

    # === Error rules catching groups of missing types by their superclasses ===

    # Generic fallback error messages for missing rules:
    def expr(self, o):
        error("Missing LaTeX formatting rule for expr type %s." % o._uflclass)

    def terminal(self, o):
        error("Missing LaTeX formatting rule for terminal type %s." % o._uflclass)

    def constant_value(self, o, component=(), derivatives=(), restriction=None):
        error("Missing LaTeX rule for constant value type %s." % o._uflclass)

    def geometric_quantity(self, o, component=(), derivatives=()):
        error("Missing LaTeX formatting rule for geometric quantity type %s." % o._uflclass)

    # Unexcepted type checks:
    def variable(self, o):
        error("Should strip away variables before formatting LaTeX code.")
        return o # or just do this if necessary

    def invalid_request(self, o, *ops):
        error("Invalid request for LaTeX formatting of a %s." % o._uflclass)
    wrapper_type = invalid_request
    index_sum = invalid_request
    indexed = invalid_request
    derivative = invalid_request
    restricted = invalid_request

    # === Formatting rules for literal constants ===

    def zero(self, o, component=(), derivatives=(), restriction=None):
        return "0" if not o.shape() else r"{\mathbf 0}"

    def int_value(self, o, component=(), derivatives=(), restriction=None):
        if derivatives:
            return self.zero(0*o)
        else:
            return "%d" % int(o)

    def float_value(self, o, component=(), derivatives=(), restriction=None):
        # Using configurable precision parameter from ufl
        if derivatives:
            return self.zero(0*o)
        else:
            return ufl.constantvalue.format_float(float(o))

    # ... The compound literals below are removed during preprocessing

    def identity(self, o):
        return r"{\mathbf I}"

    def permutation_symbol(self, o):
        return r"{\mathbf \varepsilon}"

    # === Formatting rules for geometric quantities ===

    # TODO: Add all geometric quantities here, use restriction

    def spatial_coordinate(self, o, component=(), derivatives=(), restriction=None):
        if component:
            i, = component
        else:
            i = 0
        if derivatives:
            return "x_{%d, %s}" % (i, ' '.join('%d' % d for d in derivatives))
        else:
            return "x_%d" % i

    def facet_normal(self, o, component=(), derivatives=(), restriction=None):
        if component:
            i, = component
        else:
            i = 0
        if derivatives:
            return "n_{%d, %s}" % (i, ' '.join('%d' % d for d in derivatives))
        else:
            return "n_%d" % i

    def cell_volume(self, o, component=(), derivatives=(), restriction=None):
        uflacs_assert(not component, "Expecting no component for scalar value.")
        if derivatives:
            return "0"
        else:
            return r"K_{\text{vol}}"

    def circumradius(self, o, component=(), derivatives=(), restriction=None):
        uflacs_assert(not component, "Expecting no component for scalar value.")
        if derivatives:
            return "0"
        else:
            return r"K_{\text{rad}}"

    # === Formatting rules for functions ===

    def coefficient(self, o, component=(), derivatives=(), restriction=None):
        common_name = "w"
        c = o.count()

        uflacs_assert(c >= 0, "Expecting positive count, have you preprocessed the expression?")

        name = r"\overset{%d}{%s}" % (c, common_name)

        # TODO: Use restriction

        if component:
            cstr = ' '.join('%d' % d for d in component)
        else:
            cstr = ''

        if derivatives:
            dstr = ' '.join('%d' % d for d in derivatives)
            return "%s_{%s, %s}" % (name, cstr, dstr)
        elif not component:
            return name
        else:
            return "%s_{%s}" % (name, cstr)

    def argument(self, o, component=(), derivatives=(), restriction=None):
        common_name = "v"
        c = o.number()

        name = r"\overset{%d}{%s}" % (c, common_name)

        # TODO: Use restriction

        if component:
            cstr = ' '.join('%d' % d for d in component)
        else:
            cstr = ''

        if derivatives:
            dstr = ' '.join('%d' % d for d in derivatives)
            return "%s_{%s, %s}" % (name, cstr, dstr)
        elif not component:
            return name
        else:
            return "%s_{%s}" % (name, cstr)

    # === Formatting rules for arithmetic operations ===

    def sum(self, o, *ops):
        return " + ".join(ops)

    def product(self, o, *ops):
        return " ".join(ops)

    def division(self, o, a, b):
        return r"\frac{%s}{%s}" % (a, b)

    # === Formatting rules for cmath functions ===

    def power(self, o, a, b):
        return "{%s}^{%s}" % (a, b)

    def sqrt(self, o, op):
        return "\sqrt{%s}" % (op,)

    def ln(self, o, op):
        return r"\ln(%s)" % (op,)

    def exp(self, o, op):
        return "e^{%s}" % (op,)

    def abs(self, o, op):
        return r"\|%s\|" % (op,)

    def cos(self, o, op):
        return r"\cos(%s)" % (op,)

    def sin(self, o, op):
        return r"\sin(%s)" % (op,)

    def tan(self, o, op):
        return r"\tan(%s)" % (op,)

    def cosh(self, o, op):
        return r"\cosh(%s)" % (op,)

    def sinh(self, o, op):
        return r"\sinh(%s)" % (op,)

    def tanh(self, o, op):
        return r"\tanh(%s)" % (op,)

    def acos(self, o, op):
        return r"\arccos(%s)" % (op,)

    def asin(self, o, op):
        return r"\arcsin(%s)" % (op,)

    def atan(self, o, op):
        return r"\arctan(%s)" % (op,)

    # === Formatting rules for bessel functions ===

    # TODO: Bessel functions, erf

    # === Formatting rules for conditional expressions ===

    def conditional(self, o, c, t, f):
        return r"\left{{%s} \text{if} {%s} \text{else} {%s}\right}" % (t, c, f)

    def eq(self, o, a, b):
        return r" = ".join((a, b))

    def ne(self, o, a, b):
        return r" \ne ".join((a, b))

    def le(self, o, a, b):
        return r" \le ".join((a, b))

    def ge(self, o, a, b):
        return r" \ge ".join((a, b))

    def lt(self, o, a, b):
        return r" \lt ".join((a, b))

    def gt(self, o, a, b):
        return r" \gt ".join((a, b))

    def and_condition(self, o, a, b):
        return r" \land ".join((a, b))

    def or_condition(self, o, a, b):
        return r" \lor ".join((a, b))

    def not_condition(self, o, a):
        return r" \lnot %s" % (a,)

    # === Formatting rules for restrictions ===

    def positive_restricted(self, o, a):
        return r"%s^{[+]}" % (a,) # TODO

    def negative_restricted(self, o, a):
        return r"%s^{[-]}" % (a,) # TODO


from ufl.algorithms.transformations import MultiFunction
class LatexFormatter(MultiFunction, LatexFormattingRules):
    """Default LaTeX formatter class.

    Customize by copying this and overriding rules."""
    def __init__(self):
        MultiFunction.__init__(self)

