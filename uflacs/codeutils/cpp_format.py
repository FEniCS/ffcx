
from uflacs.utils.log import info, warning, error
from uflacs.utils.assertions import uflacs_assert

import ufl
from ufl.algorithms.transformations import MultiFunction, Transformer
from ufl.algorithms import expand_derivatives, expand_indices

def build_precedence_list():
    "Builds a list of operator types by precedence order in the C language."
    # FIXME: Add all types we need here. Missing at least lifting and restriction.
    pl = []
    pl.append((ufl.classes.Conditional,))
    pl.append((ufl.classes.OrCondition,))
    pl.append((ufl.classes.AndCondition,))
    pl.append((ufl.classes.EQ, ufl.classes.NE))
    pl.append((ufl.classes.Condition,)) # <,>,<=,>=
    pl.append((ufl.classes.Sum,))
    pl.append((ufl.classes.Product, ufl.classes.Division,))
    # The highest precedence items will never need
    # parentheses around them or their operands
    pl.append((ufl.classes.Power, ufl.classes.MathFunction, ufl.classes.Abs,
               ufl.classes.Indexed, ufl.classes.SpatialDerivative,
               ufl.classes.Terminal))
    return pl

def build_precedence_map():
    from ufl.precedence import build_precedence_mapping
    pm, missing = build_precedence_mapping(build_precedence_list())
    if 0 and missing: # Enable to see which types we are missing
        print "Missing precedence levels for the types:"
        print "\n".join('  %s' % c for c in missing)
    return pm

class CodeFormatter(Transformer):
    """Language independent formatting class containing rules for
    handling indexing operators such that value and derivative
    indices are propagated to terminal handlers to be implemented
    for a particular language and target."""

    def __init__(self, language_formatter, variables):
        super(CodeFormatter, self).__init__()
        self.language_formatter = language_formatter
        self.variables = variables
        self.precedence = build_precedence_map()
        self.max_precedence = max(self.precedence.itervalues())

    def expr(self, o):
        v = self.variables.get(o)
        if v is not None:
            return v

        # Visit children and wrap in () if necessary.
        # This could be improved by considering the
        # parsing order to avoid some (), but that
        # may be language dependent? (usually left-right).
        # Keeping it simple and safe for now at least.
        ops = []
        for op in o.operands():
            opc = self.visit(op)
            po = self.precedence[o._uflclass]
            pop = self.precedence[op._uflclass]
            if po < self.max_precedence and pop <= po:
                opc = '(' + opc + ')'
            ops.append(opc)

        #ops = [self.visit(op) for op in o.operands()]
        #ops = [("(%s)" % op) for op in ops]

        return self.language_formatter(o, *ops)

    def terminal(self, o):
        v = self.variables.get(o)
        if v is not None:
            return v

        return self.language_formatter(o)

    def multi_index(self, o):
        "Expecting expand_indices to have been applied, so all indices are fixed."
        return tuple(map(int, o))

    def spatial_derivative(self, o):
        return self.spatial_derivative_component(o, ())

    def spatial_derivative_component(self, o, component):
        """Gets derivative indices and passes on control to
        either indexed or target specific terminal handler."""

        # TODO: Test variables/component/derivatives combos more!
        if 0 and self.variables:
            print ""
            print "spatial_derivative_component:"
            print self.variables
            print
            print repr(o)
            print
            print repr(self.variables.get(o))
            print
        # Use eventual given variable
        v = self.variables.get(o)
        if v is not None:
            return v

        # Note that we do not want to look for a variable for f,
        # since o represents the value of the derivative of f, not f itself.

        # o is f.dx(di)
        f, di = o.operands()

        # Sorting derivative indices, can do this because the derivatives commute
        derivatives = sorted(self.multi_index(di))

        if isinstance(f, ufl.classes.Terminal):
            # o is the derivative of a terminal expression f
            expr = f
        elif isinstance(f, ufl.classes.Indexed):
            # Since expand_indices moves Indexed in to the terminals,
            # SpatialDerivative can be outside an Indexed:
            # o is A[ci].dx(di)
            A, ci = f.operands()
            component = self.multi_index(ci)
            expr = A
        else:
            error("Invalid type %s in spatial_derivate formatter, "\
                  "have you applied expand_derivatives?" % type(o))

        # Ask the formatter to make the string
        return self.language_formatter(expr, component, derivatives)

    def indexed(self, o):
        """Gets value indices and passes on control to either
        spatial_derivative or a target specific terminal formatter."""

        # TODO: Test variables/component/derivatives combos more!
        if 0 and self.variables:
            print
            print "indexed:"
            print self.variables
            print
            print repr(o)
            print
            print repr(self.variables.get(o))
            print

        # Use eventual given variable
        v = self.variables.get(o)
        if v is not None:
            return v

        # Note that we do not want to look for a variable for
        # A, but rather for the specific component of A.
        # By only using scalar variables we keep the variables construct simple.

        # o is A[ci]
        A, ci = o.operands()
        component = self.multi_index(ci)

        if isinstance(A, ufl.classes.Terminal):
            # o is the component of a terminal A
            # Ask the formatter to make the string
            return self.language_formatter(A, component)
        elif isinstance(A, ufl.classes.SpatialDerivative):
            # A is f.dx(...)  <->  o is f.dx(...)[ci]
            # Pass on control to derivative evaluation
            return self.spatial_derivative_component(A, component)
        else:
            error("Invalid type %s in indexed formatter, "\
                  "have you applied expand_derivatives?" % type(A))


class CppLiteralFormatter(object):
    "Formatting rules for literal constants."

    def constant_value(self, o):
        error("Missing rule for constant value type %s." % o._uflclass)

    def int_value(self, o):
        return "%d" % int(o)

    def float_value(self, o):
        # Using configurable precision parameter from ufl
        return ufl.constantvalue.format_float(float(o))

    def zero(self, o):
        return "0"

class CppArithmeticFormatter(object):
    "Formatting rules for arithmetic operations."

    def sum(self, o, *ops):
        return " + ".join(ops)

    def product(self, o, *ops):
        return " * ".join(ops)

    def division(self, o, *ops):
        return " / ".join(ops)

class CppCmathFormatter(object):
    "Formatting rules for <cmath> functions."

    def power(self, o, a, b):
        return "std::pow(%s, %s)" % (a, b)

    def sqrt(self, o, op):
        return "std::sqrt(%s)" % (op,)

    def ln(self, o, op):
        return "std::log(%s)" % (op,)

    def exp(self, o, op):
        return "std::exp(%s)" % (op,)

    def abs(self, o, op):
        return "std::abs(%s)" % (op,)

    def cos(self, o, op):
        return "std::cos(%s)" % (op,)

    def sin(self, o, op):
        return "std::sin(%s)" % (op,)

    def tan(self, o, op):
        return "std::tan(%s)" % (op,)

    def acos(self, o, op):
        return "std::acos(%s)" % (op,)

    def asin(self, o, op):
        return "std::asin(%s)" % (op,)

    def atan(self, o, op):
        return "std::atan(%s)" % (op,)

class CppConditionalFormatter(object):
    "Formatting rules for conditional expressions."

    def conditional(self, o, c, t, f):
        return "%s ? %s: %s" % (c, t, f)

    def eq(self, o, a, b):
        return " == ".join((a, b))
    def ne(self, o, a, b):
        return " != ".join((a, b))
    def le(self, o, a, b):
        return " <= ".join((a, b))
    def ge(self, o, a, b):
        return " >= ".join((a, b))
    def lt(self, o, a, b):
        return " < ".join((a, b))
    def gt(self, o, a, b):
        return " > ".join((a, b))
    def and_condition(self, o, a, b):
        return " && ".join((a, b))
    def or_condition(self, o, a, b):
        return " || ".join((a, b))

class CppRestrictionFormatter(object):
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

class CppFormatterErrorRules(object):
    "Error rules catching groups of missing types by their superclasses."

    def expr(self, o):
        error("Missing C formatting rule for expr type %s." % o._uflclass)

    def terminal(self, o):
        error("Missing C formatting rule for terminal type %s." % o._uflclass)

    def variable(self, o):
        error("Should strip away variables before formatting C code.")
        return o # or just do this if necessary

    def geometric_quantity(self, o, component=(), derivatives=()):
        error("Missing C formatting rule for geometric quantity type %s." % o._uflclass)

    def algebra_operator(self, o, *ops):
        error("Missing rule for algebra operator type %s." % o._uflclass)

    def index_sum(self, o, *ops):
        error("Found an IndexSum, have you forgot to apply expand_indices?")

    def wrapper_type(self, o, *ops):
        error("Found a %s, have you forgot to apply expand_compounds or expand_indices?" % o._uflclass)

    def compound_tensor_operator(self, o, *ops):
        error("Found a %s, have you forgot to apply expand_compounds?" % o._uflclass)

    def derivative(self, o, *ops):
        error("Found a %s, have you forgot to apply expand_derivatives?" % o._uflclass)

class CppFormatterRules(MultiFunction,
                        CppFormatterErrorRules,
                        CppLiteralFormatter,
                        CppArithmeticFormatter,
                        CppCmathFormatter,
                        CppConditionalFormatter,
                        CppRestrictionFormatter):
    """Base class for target specific cpp formatter class.
    See CppTestFormatter for example of how to specialise
    for a particular target. Remember to call the constructor
    of this class from your subclass."""

    def __init__(self, target_formatter=None):
        MultiFunction.__init__(self)
        self.target_formatter = target_formatter

    # Redirect geometric quantities and form arguments to target formatter:

    def spatial_coordinate(self, o, component=(), derivatives=()):
        return self.target_formatter.spatial_coordinate(o, component, derivatives)

    def facet_normal(self, o, component=(), derivatives=()):
        return self.target_formatter.facet_normal(o, component, derivatives)

    def cell_volume(self, o, component=(), derivatives=()):
        return self.target_formatter.cell_volume(o, component, derivatives)

    def circumradius(self, o, component=(), derivatives=()):
        return self.target_formatter.circumradius(o, component, derivatives)

    def coefficient(self, o, component=(), derivatives=()):
        return self.target_formatter.coefficient(o, component, derivatives)

    def argument(self, o, component=(), derivatives=()):
        return self.target_formatter.argument(o, component, derivatives)

class CppDefaultFormatter(object):
    """Example cpp formatter class, used for the test cases.
    Override the same functions for your particular target."""
    def __init__(self):
        self.required = {}

    def require(self, o, component, derivatives, code):
        s = self.required.get(o) or {}

        key = (tuple(component), tuple(derivatives))
        oldcode = s.get(key)
        uflacs_assert((not oldcode) or (oldcode == code),
                      "Generated different code for same expression.")
        s[key] = code

        self.required[o] = s
        return code

    def spatial_coordinate(self, o, component=(), derivatives=()):
        if len(derivatives) > 1:
            return "0"

        if component:
            i, = component
        else:
            i = 0

        if derivatives:
            d, = derivatives
            return "1" if i == d else "0"
        else:
            code = "x[%d]" % i
            self.require(o, component, derivatives, code)
            return code

    def facet_normal(self, o, component=(), derivatives=()):
        if derivatives:
            return "0"

        if component:
            i, = component
        else:
            i = 0

        code = "n[%d]" % i

        self.require(o, component, derivatives, code)
        return code

    def cell_volume(self, o, component=(), derivatives=()):
        uflacs_assert(not component, "Expecting no component for scalar value.")

        if derivatives:
            return "0"

        code = "K_vol"

        self.require(o, component, derivatives, code)
        return code

    def circumradius(self, o, component=(), derivatives=()):
        uflacs_assert(not component, "Expecting no component for scalar value.")

        if derivatives:
            return "0"

        code = "K_rad"

        self.require(o, component, derivatives, code)
        return code

    def coefficient(self, o, component=(), derivatives=()):
        uflacs_assert(o.count() >= 0, "Expecting positive count, have you preprocessed the expression?")
        return self.form_argument(o, component, derivatives, 'w%d' % o.count())

    def argument(self, o, component=(), derivatives=()):
        uflacs_assert(o.count() >= 0, "Expecting positive count, have you preprocessed the expression?")
        return self.form_argument(o, component, derivatives, 'v%d' % o.count())

    def form_argument(self, o, component, derivatives, common_name):
        if derivatives:
            code = 'd%d_%s' % (len(derivatives), common_name)
            code += "".join("[%d]" % d for d in derivatives)
        else:
            code = common_name

        if component:
            code += "".join("[%d]" % c for c in component) # TODO: Or should we use a flat array?

        self.require(o, component, derivatives, code)
        return code

