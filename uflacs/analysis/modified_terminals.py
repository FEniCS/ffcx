
from six.moves import zip
from ufl.permutation import build_component_numbering
from ufl.classes import (Terminal, Grad, ReferenceGrad, Indexed, FixedIndex,
                         Restricted, PositiveRestricted, NegativeRestricted,
                         FacetAvg, CellAvg,
                         FormArgument, Coefficient, Argument, GeometricQuantity)
from ufl.sorting import sorted_expr

from ffc.log import ffc_assert, warning, error


# 
# This is THE definition of modifier types, try to use this everywhere
terminal_modifier_types = (Indexed, ReferenceGrad, Grad, Restricted, FacetAvg, CellAvg)
# 


class ModifiedTerminal(object):
    """A modified terminal expression is an object of a Terminal subtype, wrapped in terminal modifier types.

    The variables of this class are:

        expr - The original UFL expression

        terminal           - the underlying Terminal object
        global_derivatives - tuple of ints, each meaning derivative in that global direction
        local_derivatives  - tuple of ints, each meaning derivative in that local direction
        averaged           - None, 'facet' or 'cell'
        restriction        - None, '+' or '-'
        component          - tuple of ints, the global component of the Terminal
        flat_component     - single int, flattened local component of the Terminal, considering symmetry

    """
    def __init__(self, expr, terminal, global_derivatives, local_derivatives, averaged, restriction, component, flat_component):
        # The original expression
        self.expr = expr

        # The underlying terminal expression
        self.terminal = terminal

        # Components
        self.component = component
        self.flat_component = flat_component
        self.restriction = restriction

        # Derivatives
        self.global_derivatives = global_derivatives
        self.local_derivatives = local_derivatives

        # Evaluation method (alternative: { None, 'facet_midpoint', 'cell_midpoint', 'facet_avg', 'cell_avg' })
        self.averaged = averaged

    def as_tuple(self):
        t = self.terminal
        c = self.component
        gd = self.global_derivatives
        ld = self.local_derivatives
        a = self.averaged
        r = self.restriction
        return (t, c, gd, ld, a, r)

    def __hash__(self):
        return hash(self.as_tuple())

    def __eq__(self, other):
        return isinstance(other, ModifiedTerminal) and self.as_tuple() == other.as_tuple()

    def __lt__(self, other):
        return self.as_tuple() < other.as_tuple()

    def __str__(self):
        s = []
        s += ["terminal:           {0}".format(self.terminal)]
        s += ["global_derivatives: {0}".format(self.global_derivatives)]
        s += ["local_derivatives:  {0}".format(self.local_derivatives)]
        s += ["averaged:           {0}".format(self.averaged)]
        s += ["component:          {0}".format(self.component)]
        s += ["restriction:        {0}".format(self.restriction)]
        return '\n'.join(s)


def is_modified_terminal(v):
    "Check if v is a terminal or a terminal wrapped in terminal modifier types."
    while not isinstance(v, Terminal):
        if isinstance(v, terminal_modifier_types):
            v = v.operands()[0]
        else:
            return False
    return True


def strip_modified_terminal(v):
    "Extract core Terminal from a modified terminal or return None."
    while not isinstance(v, Terminal):
        if isinstance(v, terminal_modifier_types):
            v = v.operands()[0]
        else:
            return None
    return v


def analyse_modified_terminal(expr):
    """Analyse a so-called 'modified terminal' expression and return its properties in more compact form.

    A modified terminal expression is an object of a Terminal subtype, wrapped in terminal modifier types.

    The wrapper types can include 0-* Grad or ReferenceGrad objects,
    and 0-1 Restricted, 0-1 Indexed, and 0-1 FacetAvg or CellAvg objects.
    """
    # Data to determine
    component = None
    global_derivatives = []
    local_derivatives = []
    averaged = None
    restriction = None

    # Start with expr and strip away layers of modifiers
    t = expr
    while not isinstance(t, Terminal):
        if isinstance(t, Indexed):
            ffc_assert(component is None, "Got twice indexed terminal.")
            t, i = t.operands()
            ffc_assert(all(isinstance(j, FixedIndex) for j in i), "Expected only fixed indices.")
            component = [int(j) for j in i]

        elif isinstance(t, ReferenceGrad):
            ffc_assert(len(component), "Got local gradient of terminal without prior indexing.")
            local_derivatives.append(component[-1])
            component = component[:-1]
            t, = t.operands()

        elif isinstance(t, Grad):
            ffc_assert(len(component), "Got gradient of terminal without prior indexing.")
            global_derivatives.append(component[-1])
            component = component[:-1]
            t, = t.operands()

        elif isinstance(t, Restricted):
            ffc_assert(restriction is None, "Got twice restricted terminal!")
            restriction = t._side
            t, = t.operands()

        elif isinstance(t, CellAvg):
            ffc_assert(averaged is None, "Got twice averaged terminal!")
            averaged = "cell"
            t, = t.operands()

        elif isinstance(t, FacetAvg):
            ffc_assert(averaged is None, "Got twice averaged terminal!")
            averaged = "facet"
            t, = t.operands()

        elif isinstance(t, terminal_modifier_types):
            error("Missing handler for terminal modifier type %s, object is %s." % (type(t), repr(t)))

        else:
            error("Unexpected type %s object %s." % (type(t), repr(t)))

    # Make sure component is an integer tuple
    component = tuple(component) if component else ()

    # Assert that component is within the shape of the terminal (this is the global component!)
    ffc_assert(len(component) == t.rank(),
               "Length of component does not match rank of terminal.")
    ffc_assert(all(c >= 0 and c < d for c, d in zip(component, t.shape())),
               "Component indices %s are outside terminal shape %s" % (component, t.shape()))

    # Flatten component # TODO: Make the flat component local?
    symmetry = t.element().symmetry() if isinstance(t, FormArgument) else {}
    vi2si, si2vi = build_component_numbering(t.shape(), symmetry)
    flat_component = vi2si[component]
    # num_flat_components = len(si2vi)

    # Make canonical representation of derivatives
    global_derivatives = tuple(sorted(global_derivatives))
    local_derivatives = tuple(sorted(local_derivatives))

    return ModifiedTerminal(expr, t, global_derivatives, local_derivatives,
                            averaged, restriction, component, flat_component)
