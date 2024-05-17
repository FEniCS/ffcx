# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Modified terminals."""

import logging
import typing

from ufl.classes import (
    Argument,
    CellAvg,
    FacetAvg,
    FixedIndex,
    FormArgument,
    Grad,
    Indexed,
    Jacobian,
    ReferenceGrad,
    ReferenceValue,
    Restricted,
    SpatialCoordinate,
)
from ufl.permutation import build_component_numbering

logger = logging.getLogger("ffcx")


class ModifiedTerminal:
    """A modified terminal."""

    def __init__(
        self,
        expr,
        terminal,
        reference_value: bool,
        base_shape,
        base_symmetry,
        component: tuple[int, ...],
        flat_component: int,
        global_derivatives: tuple[int, ...],
        local_derivatives: tuple[int, ...],
        averaged: typing.Union[None, str],
        restriction: typing.Union[None, str],
    ):
        """Initialise.

        Args:
            expr: The original UFL expression
            terminal: the underlying Terminal object
            reference_value: whether this is represented in reference frame
            base_shape: base shape
            base_symmetry: base symmetry
            component: the global component of the Terminal
            flat_component: flattened local component of the Terminal, considering symmetry
            global_derivatives: each entry is a derivative in that global direction
            local_derivatives: each entry is a derivative in that local direction
            averaged: Entity to average over (None, 'facet' or 'cell')
            restriction: The restriction (None, '+' or '-')
        """
        # The original expression
        self.expr = expr

        # The underlying terminal expression
        self.terminal = terminal

        # Are we seeing the terminal in physical or reference frame
        self.reference_value = reference_value

        # Get the shape of the core terminal or its reference value,
        # this is the shape that component and flat_component refers to
        self.base_shape = base_shape
        self.base_symmetry = base_symmetry

        # Components
        self.component = component
        self.flat_component = flat_component

        # Derivatives
        self.global_derivatives = global_derivatives
        self.local_derivatives = local_derivatives

        # Evaluation method (alternatives: { None, 'facet_midpoint',
        #  'cell_midpoint', 'facet_avg', 'cell_avg' })
        self.averaged = averaged

        # Restriction to one cell or the other for interior facet integrals
        self.restriction = restriction

    def as_tuple(self):
        """Return a tuple with hashable values that uniquely identifies this modified terminal.

        Some of the derived variables can be omitted here as long as
        they are fully determined from the variables that are included here.
        """
        t = self.terminal  # FIXME: Terminal is not sortable...
        rv = self.reference_value
        # bs = self.base_shape
        # bsy = self.base_symmetry
        # c = self.component
        fc = self.flat_component
        gd = self.global_derivatives
        ld = self.local_derivatives
        a = self.averaged
        r = self.restriction
        return (t, rv, fc, gd, ld, a, r)

    def argument_ordering_key(self):
        """Return a key for deterministic sorting of argument vertex indices.

        The key is based on the properties of the modified terminal.
        Used in factorization but moved here for closeness with ModifiedTerminal attributes.
        """
        t = self.terminal
        assert isinstance(t, Argument)
        n = t.number()
        assert n >= 0
        p = t.part()
        rv = self.reference_value
        # bs = self.base_shape
        # bsy = self.base_symmetry
        # c = self.component
        fc = self.flat_component
        gd = self.global_derivatives
        ld = self.local_derivatives
        a = self.averaged
        r = self.restriction
        return (n, p, rv, fc, gd, ld, a, r)

    def __hash__(self):
        """Hash."""
        return hash(self.as_tuple())

    def __eq__(self, other):
        """Check equality."""
        return isinstance(other, ModifiedTerminal) and self.as_tuple() == other.as_tuple()

    def __str__(self):
        """Format as string."""
        return (
            f"terminal:           {self.terminal}\n"
            f"global_derivatives: {self.global_derivatives}\n"
            f"local_derivatives:  {self.local_derivatives}\n"
            f"averaged:           {self.averaged}\n"
            f"component:          {self.component}\n"
            f"restriction:        {self.restriction}"
        )


def is_modified_terminal(v):
    """Check if v is a terminal or a terminal wrapped in terminal modifier types."""
    while not v._ufl_is_terminal_:
        if v._ufl_is_terminal_modifier_:
            v = v.ufl_operands[0]
        else:
            return False
    return True


def strip_modified_terminal(v):
    """Extract core Terminal from a modified terminal or return None."""
    while not v._ufl_is_terminal_:
        if v._ufl_is_terminal_modifier_:
            v = v.ufl_operands[0]
        else:
            return None
    return v


def analyse_modified_terminal(expr):
    """Analyse a so-called 'modified terminal' expression.

    Return its properties in more compact form as a ModifiedTerminal object.

    A modified terminal expression is an object of a Terminal subtype,
    wrapped in terminal modifier types.

    The wrapper types can include 0-* Grad or ReferenceGrad objects,
    and 0-1 ReferenceValue, 0-1 Restricted, 0-1 Indexed,
    and 0-1 FacetAvg or CellAvg objects.
    """
    # Data to determine
    component = None
    global_derivatives = []
    local_derivatives = []
    reference_value = None
    restriction = None
    averaged = None

    # Start with expr and strip away layers of modifiers
    t = expr
    while not t._ufl_is_terminal_:
        if isinstance(t, Indexed):
            if component is not None:
                raise RuntimeError("Got twice indexed terminal.")

            t, i = t.ufl_operands
            component = [int(j) for j in i]

            if not all(isinstance(j, FixedIndex) for j in i):
                raise RuntimeError("Expected only fixed indices.")

        elif isinstance(t, ReferenceValue):
            if reference_value is not None:
                raise RuntimeError("Got twice pulled back terminal!")

            (t,) = t.ufl_operands
            reference_value = True

        elif isinstance(t, ReferenceGrad):
            if not component:  # covers None or ()
                raise RuntimeError("Got local gradient of terminal without prior indexing.")

            (t,) = t.ufl_operands
            local_derivatives.append(component[-1])
            component = component[:-1]

        elif isinstance(t, Grad):
            if not component:  # covers None or ()
                raise RuntimeError("Got local gradient of terminal without prior indexing.")

            (t,) = t.ufl_operands
            global_derivatives.append(component[-1])
            component = component[:-1]

        elif isinstance(t, Restricted):
            if restriction is not None:
                raise RuntimeError("Got twice restricted terminal!")

            restriction = t._side
            (t,) = t.ufl_operands

        elif isinstance(t, CellAvg):
            if averaged is not None:
                raise RuntimeError("Got twice averaged terminal!")

            (t,) = t.ufl_operands
            averaged = "cell"

        elif isinstance(t, FacetAvg):
            if averaged is not None:
                raise RuntimeError("Got twice averaged terminal!")

            (t,) = t.ufl_operands
            averaged = "facet"

        elif t._ufl_terminal_modifiers_:
            raise RuntimeError(
                f"Missing handler for terminal modifier type {type(t)}, object is {t!r}."
            )
        else:
            raise RuntimeError(f"Unexpected type {type(t)} object {t}.")

    # Make canonical representation of derivatives
    global_derivatives = tuple(sorted(global_derivatives))
    local_derivatives = tuple(sorted(local_derivatives))

    # Make reference_value true or false
    reference_value = reference_value or False

    # Consistency check
    if isinstance(t, (SpatialCoordinate, Jacobian)):
        pass
    else:
        if local_derivatives and not reference_value:
            raise RuntimeError("Local derivatives of non-local value is not legal.")
        if global_derivatives and reference_value:
            raise RuntimeError("Global derivatives of local value is not legal.")

    # Make sure component is an integer tuple
    if component is None:
        component = ()
    else:
        component = tuple(component)

    # Get the shape of the core terminal or its reference value, this is
    # the shape that component refers to
    if isinstance(t, FormArgument):
        element = t.ufl_function_space().ufl_element()
        if reference_value:
            # Ignoring symmetry, assuming already applied in conversion
            # to reference frame
            base_symmetry = {}
            base_shape = element.reference_value_shape
        else:
            base_symmetry = element.symmetry()
            base_shape = t.ufl_shape
    else:
        base_symmetry = {}
        base_shape = t.ufl_shape

    # Assert that component is within the shape of the (reference)
    # terminal
    if len(component) != len(base_shape):
        raise RuntimeError("Length of component does not match rank of (reference) terminal.")
    if not all(c >= 0 and c < d for c, d in zip(component, base_shape)):
        raise RuntimeError("Component indices {component} are outside value shape {base_shape}.")

    # Flatten component
    vi2si, _ = build_component_numbering(base_shape, base_symmetry)
    flat_component = vi2si[component]

    return ModifiedTerminal(
        expr,
        t,
        reference_value,
        base_shape,
        base_symmetry,
        component,
        flat_component,
        global_derivatives,
        local_derivatives,
        averaged,
        restriction,
    )
