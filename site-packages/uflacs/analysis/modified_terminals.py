
from ufl.classes import (Terminal, Grad, Indexed, FixedIndex,
                         Restricted, PositiveRestricted, NegativeRestricted,
                         FacetAvg, CellAvg,
                         Coefficient, Argument, GeometricQuantity)
from ufl.sorting import sorted_expr

from uflacs.utils.log import uflacs_assert, warning, error

# TODO: Add FacetAvg and CellAvg to modifiers everywhere relevant and handle in table extraction
# TODO: Make this more robust by looping like analyse_modified_terminal, currently assumes that transformations have been applied.
def is_modified_terminal(v):
    _accepted_types = (Terminal, Grad, Restricted, FacetAvg, CellAvg)
    return (isinstance(v, _accepted_types)
            or (isinstance(v, Indexed) and isinstance(v.operands()[0], _accepted_types)))

class ModifiedTerminal(object):
    def __init__(self, terminal, derivatives, averaged, restriction, component):
        self.terminal = terminal
        self.derivatives = derivatives
        self.averaged = averaged
        self.restriction = restriction
        self.component = component

terminal_modifier_types = (Grad, Restricted, Indexed, FacetAvg, CellAvg)
def analyse_modified_terminal2(o, form_argument_mapping={}):
    """Analyse a so-called 'modified terminal' expression and return its properties in more compact form.

    A modified terminal expression is:
    an object of a Terminal subtype,
    wrapped in 0-* Grad objects,
    wrapped in 0-1 Restricted object,
    wrapped in 0-1 Indexed object.

    The returned values are:

    (terminal, component, derivatives, restriction)

    # TODO: Explain the format of these values for future reference.
    """
    t = o
    component = None
    derivatives = []
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

        elif isinstance(t, Grad):
            uflacs_assert(len(component), "Got gradient of terminal without prior indexing.")
            derivatives.append(component[-1])
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

    t = form_argument_mapping.get(t,t)
    component = tuple(component) if component else ()
    derivatives = tuple(sorted(derivatives))

    uflacs_assert(len(component) == t.rank(),
                  "Length of component does not match rank of terminal.")
    uflacs_assert(all(c >= 0 and c < d for c,d in zip(component, t.shape())),
                  "Component indices %s are outside terminal shape %s" % (component, t.shape()))

    # FIXME: Return mt and update all callers, then use averaged state
    mt = ModifiedTerminal(t, derivatives, averaged, restriction, component)
    return mt

def analyse_modified_terminal(o, form_argument_mapping={}): # FIXME: Temporary wrapper to transition from tuple to mt struct
    mt = analyse_modified_terminal2(o, form_argument_mapping)
    return (mt.terminal, mt.component, mt.derivatives, mt.restriction)
