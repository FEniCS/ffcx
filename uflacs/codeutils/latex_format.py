
from uflacs.utils.log import info, warning, error
from uflacs.utils.assertions import uflacs_assert

import ufl
from ufl.algorithms.transformations import MultiFunction

class LatexLiteralFormatter(object):
    "Formatting rules for literal constants."

    def constant_value(self, o):
        error("Missing rule for constant value type %s." % o._uflclass)

    def int_value(self, o):
        return "%d" % int(o)

    def float_value(self, o):
        # Using configurable precision parameter from ufl
        return ufl.constantvalue.format_float(float(o))

    def zero(self, o):
        return "0" if not o.shape() else r"{\mathbf 0}"

    def identity(self, o):
        return r"{\mathbf I}"

    def permutation_symbol(self, o):
        return r"{\mathbf \varepsilon}"

class LatexGeometryFormatter(object):

    def spatial_coordinate(self, o, component=(), derivatives=()):
        if component:
            i, = component
        else:
            i = 0
        if derivatives:
            return "x_{%d, %s}" % (i, ' '.join('%d' % d for d in derivatives))
        else:
            return "x_%d" % i

    def facet_normal(self, o, component=(), derivatives=()):
        if component:
            i, = component
        else:
            i = 0
        if derivatives:
            return "n_{%d, %s}" % (i, ' '.join('%d' % d for d in derivatives))
        else:
            return "n_%d" % i

    def cell_volume(self, o, component=(), derivatives=()):
        uflacs_assert(not component, "Expecting no component for scalar value.")
        if derivatives:
            return "0"
        else:
            return r"K_{\text{vol}}"

    def circumradius(self, o, component=(), derivatives=()):
        uflacs_assert(not component, "Expecting no component for scalar value.")
        if derivatives:
            return "0"
        else:
            return r"K_{\text{rad}}"

class LatexFormArgumentFormatter(object):

    def coefficient(self, o, component=(), derivatives=()):
        uflacs_assert(o.count() >= 0, "Expecting positive count, have you preprocessed the expression?")
        return self.form_argument(o, component, derivatives, 'w')

    def argument(self, o, component=(), derivatives=()):
        uflacs_assert(o.count() >= 0, "Expecting positive count, have you preprocessed the expression?")
        return self.form_argument(o, component, derivatives, 'v')

    def form_argument(self, o, component, derivatives, common_name):
        name = r"\overset{%d}{%s}" % (o.count(), common_name)
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

class LatexArithmeticFormatter(object):
    "Formatting rules for arithmetic operations."

    def sum(self, o, *ops):
        return " + ".join(ops)

    def product(self, o, *ops):
        return " ".join(ops)

    def division(self, o, a, b):
        return r"\frac{%s}{%s}" % (a, b)

class LatexCmathFormatter(object):
    "Formatting rules for <cmath> functions."

    def power(self, o, a, b):
        return "pow(%s, %s)" % (a, b) # TODO
        #return "{%s}^{%s}" % (a, b)

    def sqrt(self, o, op):
        return "sqrt(%s)" % (op,) # TODO
        #return r"%s^{\frac 1 2}" % (op,)

    def ln(self, o, op):
        return r"\ln(%s)" % (op,)

    def exp(self, o, op):
        #return "e^{%s}" % (op,)
        return "exp(%s)" % (op,) # TODO

    def abs(self, o, op):
        return r"\|%s\|" % (op,)

    def cos(self, o, op):
        return r"\cos(%s)" % (op,)

    def sin(self, o, op):
        return r"\sin(%s)" % (op,)

    def tan(self, o, op):
        return r"\tan(%s)" % (op,)

    def acos(self, o, op):
        return r"\arccos(%s)" % (op,)

    def asin(self, o, op):
        return r"\arcsin(%s)" % (op,)

    def atan(self, o, op):
        return r"\arctan(%s)" % (op,)

class LatexConditionalFormatter(object):
    "Formatting rules for conditional expressions."

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

class LatexRestrictionFormatter(object):
    "Not sure how to handle this or wether it makes sense to have a formatter for it."
    def restricted(self, o, *ops):
        raise NotImplemented
    def positive_restricted(self, o, *ops):
        raise NotImplemented
    def negative_restricted(self, o, *ops):
        raise NotImplemented

    def lifting_result(self, o, *ops):
        raise NotImplemented
    def lifting_operator_result(self, o, *ops):
        raise NotImplemented
    def lifting_function_result(self, o, *ops):
        raise NotImplemented
    def lifting_operator(self, o, *ops):
        raise NotImplemented
    def lifting_function(self, o, *ops):
        raise NotImplemented

class LatexFormatterErrorRules(object):
    "Error rules catching groups of missing types by their superclasses."

    def expr(self, o):
        error("Missing LaTeX formatting rule for expr type %s." % o._uflclass)

    def terminal(self, o):
        error("Missing LaTeX formatting rule for terminal type %s." % o._uflclass)

    def variable(self, o):
        error("Should strip away variables before formatting LaTeX code.")
        return o # or just do this if necessary

    def geometric_quantity(self, o, component=(), derivatives=()):
        error("Missing LaTeX formatting rule for geometric quantity type %s." % o._uflclass)

    def algebra_operator(self, o, *ops):
        error("Missing LaTeX formatting rule for algebra operator type %s." % o._uflclass)

    def index_sum(self, o, *ops):
        error("Found an IndexSum, have you forgot to apply expand_indices?")

    def wrapper_type(self, o, *ops):
        error("Found a %s, have you forgot to apply expand_compounds or expand_indices?" % o._uflclass)

    def compound_tensor_operator(self, o, *ops):
        error("Found a %s, have you forgot to apply expand_compounds?" % o._uflclass)

    def derivative(self, o, *ops):
        error("Found a %s, have you forgot to apply expand_derivatives?" % o._uflclass)

class LatexFormatterRules(MultiFunction,
                          LatexFormatterErrorRules,
                          LatexLiteralFormatter,
                          LatexGeometryFormatter,
                          LatexFormArgumentFormatter,
                          LatexArithmeticFormatter,
                          LatexCmathFormatter,
                          LatexConditionalFormatter,
                          LatexRestrictionFormatter):
    """LaTeX formatter class."""
    def __init__(self):
        MultiFunction.__init__(self)
