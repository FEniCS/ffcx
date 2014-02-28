
from ufl.classes import (Terminal, Grad, LocalGrad, Indexed, FixedIndex,
                         Restricted, PositiveRestricted, NegativeRestricted,
                         FacetAvg, CellAvg,
                         Coefficient, Argument, GeometricQuantity)
from ufl.sorting import sorted_expr

from uflacs.utils.log import uflacs_assert, warning, error


# TODO: Add LocalGrad, FacetAvg and CellAvg to modifiers everywhere relevant and handle in table extraction


# This is THE definition of modifier types, try to use this everywhere
terminal_modifier_types = (LocalGrad, Grad, Restricted, Indexed, FacetAvg, CellAvg)

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


class ModifiedTerminal(object):
    """A modified terminal expression is an object of a Terminal subtype, wrapped in terminal modifier types.

    The variables of this class are:

        terminal - the Terminal object
        global_derivatives - TODO: Explain
        local_derivatives  - TODO: Explain
        averages    - TODO: Explain
        restriction - TODO: Explain
        component   - TODO: Explain
    """
    def __init__(self, terminal, global_derivatives, local_derivatives, averaged, restriction, component):
        self.terminal = terminal
        self.global_derivatives = global_derivatives
        self.local_derivatives = local_derivatives
        self.averaged = averaged
        self.restriction = restriction
        self.component = component


def analyse_modified_terminal2(o):
    """Analyse a so-called 'modified terminal' expression and return its properties in more compact form.

    A modified terminal expression is an object of a Terminal subtype, wrapped in terminal modifier types.

    The wrapper types can include 0-* Grad or LocalGrad objects,
    and 0-1 Restricted, 0-1 Indexed, and 0-1 FacetAvg or CellAvg objects.
    """
    t = o
    component = None
    global_derivatives = []
    local_derivatives = []
    restriction = None
    averaged = None
    while not isinstance(t, Terminal):
        if not isinstance(t, terminal_modifier_types):
            error("Unexpected type %s object %s." % (type(t), repr(t)))

        if isinstance(t, Indexed):
            uflacs_assert(component is None, "Got twice indexed terminal.")
            t, i = t.operands()
            uflacs_assert(all(isinstance(j, FixedIndex) for j in i), "Expected only fixed indices.")
            component = [int(j) for j in i]

        elif isinstance(t, LocalGrad):
            uflacs_assert(len(component), "Got local gradient of terminal without prior indexing.")
            local_derivatives.append(component[-1])
            component = component[:-1]
            t, = t.operands()

        elif isinstance(t, Grad):
            uflacs_assert(len(component), "Got gradient of terminal without prior indexing.")
            global_derivatives.append(component[-1])
            component = component[:-1]
            t, = t.operands()

        elif isinstance(t, Restricted):
            uflacs_assert(restriction is None, "Got twice restricted terminal!")
            restriction = t._side
            t, = t.operands()

        elif isinstance(t, CellAvg):
            uflacs_assert(averaged is None, "Got twice averaged terminal!")
            averaged = "cell"
            t, = t.operands()

        elif isinstance(t, FacetAvg):
            uflacs_assert(averaged is None, "Got twice averaged terminal!")
            averaged = "facet"
            t, = t.operands()

    component = tuple(component) if component else ()
    global_derivatives = tuple(sorted(global_derivatives))
    local_derivatives = tuple(sorted(local_derivatives))

    uflacs_assert(len(component) == t.rank(),
                  "Length of component does not match rank of terminal.")
    uflacs_assert(all(c >= 0 and c < d for c,d in zip(component, t.shape())),
                  "Component indices %s are outside terminal shape %s" % (component, t.shape()))

    return ModifiedTerminal(t, global_derivatives, local_derivatives, averaged, restriction, component)


def analyse_modified_terminal(o): # FIXME: Temporary wrapper to transition from tuple to mt struct
    mt = analyse_modified_terminal2(o)
    return (mt.terminal, mt.component, mt.derivatives, mt.restriction)

