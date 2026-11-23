# Copyright (C) 2011-2026 Martin Sandve Alnæs, Jørgen S. Dokken and Paul T. Kühner
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Modified terminals."""

import logging
from typing import TypeAlias

from ufl.algorithms.analysis import extract_type
from ufl.classes import (
    Argument,
    CellAvg,
    Coefficient,
    Expr,
    FacetAvg,
    FixedIndex,
    Grad,
    Indexed,
    Interpolate,
    Jacobian,
    ReferenceGrad,
    ReferenceValue,
    Restricted,
    SpatialCoordinate,
    Terminal,
)
from ufl.utils.indexflattening import flatten_multiindex, shape_to_strides

logger = logging.getLogger("ffcx")

CoreTerminal: TypeAlias = Terminal | Interpolate
"""The core of a modified terminal.
A :class:`ufl.core.terminal.Terminal`, or an :class:`ufl.Interpolate`.
"""


class ModifiedTerminal:
    """A modified terminal."""

    def __init__(
        self,
        expr: Expr,
        terminal: CoreTerminal,
        reference_value: bool,
        base_shape: tuple[int, ...],
        component: tuple[int, ...],
        flat_component: int,
        global_derivatives: tuple[int, ...],
        local_derivatives: tuple[int, ...],
        averaged: str | None,
        restriction: str | None,
    ):
        """Initialise.

        Args:
            expr: The original :class:`ufl.core.expr.Expr`
            terminal: the core terminal, a :class:`ufl.core.terminal.Terminal` or an
                :class:`ufl.Interpolate`
            reference_value: whether this is represented in reference frame
            base_shape: base shape
            component: the global component of the Terminal
            flat_component: flattened local component of the Terminal
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
        arguments = modified_terminal_arguments(self.terminal)
        assert arguments
        # An interpolation may stand for several arguments; order by the first.
        n = arguments[0].number()
        assert n >= 0
        p = arguments[0].part()
        rv = self.reference_value
        # bs = self.base_shape
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


def modified_terminal_arguments(terminal: CoreTerminal) -> tuple[Argument, ...]:
    """Get the arguments that the core terminal of a modified terminal depends on.

    Each argument of a form owns one axis of the local tensor, ordered by argument
    number. The axis runs over the argument's basis functions on the cell of the
    integration entity, or on both cells for an interior facet integral.

    An argument depends on itself and contributes along its own axis.
    An :class:`ufl.Interpolate` depends on the arguments of the expression it
    interpolates and contributes along their axes, without adding an axis of its own.
    The :class:`ufl.Coargument` in its dual slot is contracted by the surrounding
    expression, so it is not among them. Any other terminal depends on none.

    Args:
        terminal: The core terminal of a modified terminal.

    Returns:
        The arguments, ordered by number and part.
    """
    arguments = extract_type(terminal, Argument, base_form_op_as_expr=True)
    return tuple(
        sorted(arguments, key=lambda a: (a.number(), -1 if a.part() is None else a.part()))
    )


def is_modified_terminal(v: Expr) -> bool:
    """Check if v is a terminal or a terminal wrapped in terminal modifier types.

    An interpolation counts as a terminal: the expression it interpolates is
    evaluated separately, at the interpolation points of its target element.
    """
    while not v._ufl_is_terminal_:
        if isinstance(v, Interpolate):
            return True
        elif v._ufl_is_terminal_modifier_:
            v = v.ufl_operands[0]
        else:
            return False
    return True


def strip_modified_terminal(v: Expr) -> CoreTerminal | None:
    """Extract the core terminal from a modified terminal, or return None."""
    while not v._ufl_is_terminal_:
        if isinstance(v, Interpolate):
            return v
        elif v._ufl_is_terminal_modifier_:
            v = v.ufl_operands[0]
        else:
            return None
    assert isinstance(v, Terminal)
    return v


def analyse_modified_terminal(expr: Expr) -> ModifiedTerminal:
    """Analyse a so-called 'modified terminal' expression.

    Return its properties in more compact form as a ModifiedTerminal object.

    A modified terminal expression is a core terminal, a
    :class:`ufl.core.terminal.Terminal` or an :class:`ufl.Interpolate`, wrapped in
    terminal modifier types.

    The wrapper types can include 0-* Grad or ReferenceGrad objects,
    and 0-1 ReferenceValue, 0-1 Restricted, 0-1 Indexed,
    and 0-1 FacetAvg or CellAvg objects.
    """
    # Data to determine
    component: list[int] | None = None
    global_derivatives: list[int] = []
    local_derivatives: list[int] = []
    reference_value = None
    restriction = None
    averaged = None

    # Start with expr and strip away layers of modifiers
    t = expr
    while not (t._ufl_is_terminal_ or isinstance(t, Interpolate)):
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
    assert isinstance(t, Terminal | Interpolate)

    # Make canonical representation of derivatives
    sorted_global_derivatives = tuple(sorted(global_derivatives))
    sorted_local_derivatives = tuple(sorted(local_derivatives))

    # Make reference_value true or false
    reference_value = reference_value or False

    # Consistency check
    if isinstance(t, SpatialCoordinate | Jacobian):
        pass
    else:
        if local_derivatives and not reference_value:
            raise RuntimeError("Local derivatives of non-local value is not legal.")
        if global_derivatives and reference_value:
            raise RuntimeError("Global derivatives of local value is not legal.")

    # Make sure component is an integer tuple
    component_tuple: tuple[int, ...] = () if component is None else tuple(component)

    # Get the shape of the core terminal or its reference value, this is
    # the shape that component refers to. Form arguments are read in
    # reference value, with the pull backs of their elements applied and
    # any symmetry with them, so their components carry no symmetry.
    if isinstance(t, Argument | Coefficient | Interpolate):
        if not reference_value:
            raise RuntimeError(f"Expected the reference value of {t}: pull backs not applied.")
        # An interpolation is a finite element field of its target space.
        base_shape = t.ufl_element().reference_value_shape
    else:
        base_shape = t.ufl_shape

    # Assert that component is within the shape of the (reference)
    # terminal
    if len(component_tuple) != len(base_shape):
        raise RuntimeError("Length of component does not match rank of (reference) terminal.")
    if not all(c >= 0 and c < d for c, d in zip(component_tuple, base_shape)):
        raise RuntimeError(
            f"Component indices {component_tuple} are outside value shape {base_shape}."
        )

    # Flatten component
    flat_component = flatten_multiindex(component_tuple, shape_to_strides(base_shape))

    return ModifiedTerminal(
        expr,
        t,
        reference_value,
        base_shape,
        component_tuple,
        flat_component,
        sorted_global_derivatives,
        sorted_local_derivatives,
        averaged,
        restriction,
    )
