# Copyright (C) 2007-2026 Anders Logg, Martin Alnaes, Kristian B. Oelgaard,
#                         Michal Habera and others
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Compiler stage 1: Analysis.

This module implements the analysis/preprocessing of variational forms,
including automatic selection of elements, degrees and form
representation type.
"""

from __future__ import annotations

import logging
import typing
from functools import singledispatchmethod

if typing.TYPE_CHECKING:
    from ufl.algorithms.formdata import FormData

import basix.ufl
import numpy as np
import numpy.typing as npt
import ufl.algorithms
from ufl.algorithms.apply_algebra_lowering import apply_algebra_lowering
from ufl.algorithms.apply_derivatives import apply_coordinate_derivatives, apply_derivatives
from ufl.algorithms.apply_function_pullbacks import (
    apply_function_pullbacks,
    apply_interpolate_pullbacks,
)
from ufl.algorithms.apply_geometry_lowering import apply_geometry_lowering
from ufl.algorithms.apply_integral_scaling import apply_integral_scaling
from ufl.algorithms.compute_form_data import attach_estimated_degrees, preprocess_form

# See TODOs at the call sites of these below:
from ufl.algorithms.domain_analysis import (
    build_integral_data,
    group_form_integrals,
)
from ufl.algorithms.formdata import FormData
from ufl.algorithms.remove_complex_nodes import remove_complex_nodes
from ufl.algorithms.remove_component_tensors import remove_component_tensors
from ufl.corealg.dag_traverser import DAGTraverser

logger = logging.getLogger("ffcx")


def apply_push_forward(
    expr: ufl.core.expr.Expr, domain: ufl.Mesh, pullback: ufl.pullback.AbstractPullback
) -> ufl.core.expr.Expr:
    """Apply push forward to an expression."""
    J = ufl.Jacobian(domain)
    J_T = J.T
    detJ = ufl.JacobianDeterminant(domain)
    invJ = ufl.JacobianInverse(domain)
    invJ_T = invJ.T
    pushed_forward_expr: ufl.core.expr.Expr
    match pullback:
        case ufl.pullback.L2Piola():
            pushed_forward_expr = detJ * expr
        case ufl.pullback.CovariantPiola():
            pushed_forward_expr = ufl.dot(J_T, expr)
        case ufl.pullback.ContravariantPiola():
            pushed_forward_expr = detJ * ufl.dot(invJ, expr)
        case ufl.pullback.DoubleCovariantPiola():
            pushed_forward_expr = ufl.dot(J_T, ufl.dot(expr, J))
        case ufl.pullback.DoubleContravariantPiola():
            pushed_forward_expr = detJ * ufl.dot(invJ, ufl.dot(expr, invJ_T))
        case ufl.pullback.CovariantContravariantPiola():
            pushed_forward_expr = detJ * ufl.dot(J_T, ufl.dot(expr, invJ_T))
        case ufl.pullback.IdentityPullback():
            pushed_forward_expr = expr
        case _:
            raise NotImplementedError(f"Pullback {pullback} not supported.")
    return pushed_forward_expr


class UFLData(typing.NamedTuple):
    """UFL data."""

    #: FormData objects
    form_data: tuple[FormData, ...]
    #: Unique elements across all forms and expressions
    unique_elements: list[basix.ufl._ElementBase]
    #: Mapping to unique numbers for all elements
    element_numbers: dict[basix.ufl._ElementBase, int]
    #: Unique coordinate elements across all forms and expressions
    unique_coordinate_elements: list[basix.ufl._ElementBase]
    #: List of all expressions after post-processing, with evaluation points and original expression
    expressions: list[tuple[ufl.core.expr.Expr, npt.NDArray[np.floating], ufl.core.expr.Expr]]


def _interpolations(data: FormData) -> list[ProxyCoefficient | ufl.Interpolate]:
    """Collect what needs an expression kernel evaluated at interpolation points.

    A proxy coefficient always needs one: its degrees of freedom are built per
    cell by evaluating the interpolated expression. An interpolation that
    survived into the integrands, which is one holding an argument, needs one
    too, unless its expression is the argument itself and the table is known at
    compile time.

    Args:
        data: Form data to collect from.
    """
    collected: list[ProxyCoefficient | ufl.Interpolate] = [
        coefficient
        for coefficient in data.reduced_coefficients
        if isinstance(coefficient, ProxyCoefficient)
    ]
    seen = set()
    for integral_data in data.integral_data:
        for integral in integral_data.integrals:
            for node in ufl.corealg.traversal.unique_pre_traversal(integral.integrand()):
                if (
                    isinstance(node, ufl.Interpolate)
                    and not isinstance(node.ufl_operands[0], ufl.Argument)
                    and node not in seen
                ):
                    seen.add(node)
                    collected.append(node)
    return collected


def interpolated_expression(
    interpolation: ProxyCoefficient | ufl.Interpolate,
) -> ufl.core.expr.Expr:
    """Get the expression that an interpolation or proxy coefficient evaluates.

    Args:
        interpolation: A `ProxyCoefficient` or a `ufl.Interpolate`.
    """
    if isinstance(interpolation, ProxyCoefficient):
        return interpolation.operand
    (operand,) = interpolation.ufl_operands
    return operand


def interpolation_arguments(
    interpolation: ProxyCoefficient | ufl.Interpolate,
) -> tuple[ufl.Argument, ...]:
    """Get the arguments that an interpolated expression is linear in, by number.

    Args:
        interpolation: A `ProxyCoefficient` or a `ufl.Interpolate`.
    """
    arguments = ufl.algorithms.extract_arguments(interpolated_expression(interpolation))
    return tuple(sorted(arguments, key=lambda argument: argument.number()))


def interpolation_dof_elements(
    interpolation: ufl.Interpolate,
) -> tuple[basix.ufl._ElementBase, ...]:
    """Get the elements whose degrees of freedom an interpolation's table is indexed by.

    Ordered by argument number, so there are two of them for the interpolation
    of an expression bilinear in two arguments.

    Args:
        interpolation: A `ufl.Interpolate`.
    """
    return tuple(
        argument.ufl_function_space().ufl_element()
        for argument in interpolation_arguments(interpolation)
    )


def interpolated_argument(interpolation: ufl.Interpolate) -> ufl.Argument | None:
    """Get the argument that an interpolation maps directly, if it maps one.

    `apply_interpolate_pullbacks` maps the operand onto the reference cell, so
    interpolating an argument itself leaves the argument under the target
    element's inverse pull back rather than on its own.

    Args:
        interpolation: A `ufl.Interpolate`.

    Returns:
        The argument, or None if the operand is not one mapped that way.
    """
    operand = interpolated_expression(interpolation)
    arguments = interpolation_arguments(interpolation)
    if len(arguments) != 1:
        return None
    (argument,) = arguments
    space = interpolation_target_space(interpolation)
    target = space.ufl_element()
    source = argument.ufl_function_space().ufl_element()
    domain = ufl.domain.extract_unique_domain(argument) or space.ufl_domain()
    # The argument as it stands in the operand. `apply_interpolate_pullbacks`
    # maps the operand onto the reference cell of the target element, and
    # `apply_function_pullbacks` later represents the argument itself in
    # reference value, through its own pull back. Both forms are met: the
    # interpolations are checked between the two passes, and classified after.
    mapped = (
        target.pullback.apply_inverse(argument, domain),
        target.pullback.apply_inverse(
            source.pullback.apply(ufl.classes.ReferenceValue(argument), domain), domain
        ),
    )
    return argument if any(_same_expression(operand, m) for m in mapped) else None


def _same_expression(a: ufl.core.expr.Expr, b: ufl.core.expr.Expr) -> bool:
    """Compare two expressions up to lowering and index numbering.

    The operand of an interpolation is compared against one built here, which
    has not been through the passes the form has: the geometry of a pull back is
    lowered in the form, and each pull back numbers its free indices from a
    global counter.

    Args:
        a: An expression.
        b: An expression to compare it against.
    """

    def normalise(expression: ufl.core.expr.Expr) -> ufl.core.expr.Expr:
        preserve = (ufl.classes.Jacobian,)
        for _ in range(2):
            expression = apply_geometry_lowering(expression, preserve)
            expression = apply_derivatives(expression)
        return ufl.algorithms.renumbering.renumber_indices(expression)

    return bool(normalise(a) == normalise(b))


def interpolation_has_runtime_table(interpolation: ufl.Interpolate) -> bool:
    """Whether an interpolation's element table has to be built for each cell.

    Interpolating an argument itself is a fixed map between the reference
    elements, so the table is known at compile time. Interpolating an expression
    that is merely linear in the argument, as differentiating a non-linear
    interpolated expression produces, brings in coefficients, so the table
    depends on the cell.

    Args:
        interpolation: A `ufl.Interpolate`.
    """
    return interpolated_argument(interpolation) is None


def check_interpolation(interpolation: ufl.Interpolate) -> None:
    """Check that an interpolation holding an argument can be turned into a table.

    Args:
        interpolation: A `ufl.Interpolate`.
    """
    target = interpolation_target_space(interpolation).ufl_element()
    argument = interpolated_argument(interpolation)
    if argument is not None:
        # The argument's reference basis is pushed straight through the
        # interpolation, so the table is known at compile time. That needs the
        # Jacobian factors of the two pull backs to cancel and the block
        # structures to line up.
        source = argument.ufl_function_space().ufl_element()
        if source.pullback != target.pullback:
            raise NotImplementedError(
                f"Interpolation of an argument from {source.pullback} to "
                f"{target.pullback} is not supported."
            )
        if source.block_size != target.block_size:
            raise NotImplementedError(
                f"Interpolation of an argument from a block size {source.block_size} "
                f"space into a block size {target.block_size} space is not supported."
            )
    else:
        # The table has to be built for each cell from an expression kernel, see
        # `ffcx.codegeneration.integral_generator.IntegralGenerator`.
        if target.pullback != ufl.pullback.identity_pullback:
            raise NotImplementedError(
                f"Interpolating an expression into a {target.pullback} space is "
                "only supported when the expression is the argument itself."
            )
        if target.block_size != target.reference_value_size:
            raise NotImplementedError(
                "Interpolating an expression into a non-blocked vector valued "
                "space is only supported when the expression is the argument itself."
            )


def interpolation_target_space(
    interpolation: ProxyCoefficient | ufl.Interpolate,
) -> ufl.FunctionSpace:
    """Get the function space that an interpolation lands in.

    Args:
        interpolation: A `ProxyCoefficient` or a `ufl.Interpolate`.
    """
    if isinstance(interpolation, ProxyCoefficient):
        return interpolation.ufl_function_space()
    return interpolation.target_space()


def analyze_ufl_objects(
    ufl_objects: list[
        ufl.form.Form
        | basix.ufl._ElementBase
        | ufl.Mesh
        | tuple[ufl.core.expr.Expr, npt.NDArray[np.floating]]
    ],
    scalar_type: npt.DTypeLike,
) -> UFLData:
    """Analyze ufl object(s).

    Args:
        ufl_objects: UFL objects
        scalar_type: Scalar type that should be used for the analysis

    Returns:
        A named tuple :class:`UFLData`.
    """
    logger.info(79 * "*")
    logger.info("Compiler stage 1: Analyzing UFL objects")
    logger.info(79 * "*")

    elements: list[basix.ufl._ElementBase] = []
    coordinate_elements: list[basix.ufl._ElementBase] = []

    # Group objects by types
    forms: list[ufl.form.Form] = []
    expressions: list[tuple[ufl.core.expr.Expr, npt.NDArray[np.floating]]] = []
    processed_expressions: list[
        tuple[ufl.core.expr.Expr, npt.NDArray[np.floating], ufl.core.expr.Expr]
    ] = []

    for ufl_object in ufl_objects:
        if isinstance(ufl_object, ufl.form.Form):
            forms.append(ufl_object)
        elif isinstance(ufl_object, ufl.AbstractFiniteElement):
            elements.append(ufl_object)
        elif isinstance(ufl_object, ufl.Mesh):
            coordinate_elements.append(ufl_object.ufl_coordinate_element())
        elif isinstance(ufl_object[0], ufl.core.expr.Expr):
            original_expression = ufl_object[0]
            points = np.asarray(ufl_object[1])
            expressions.append((original_expression, points))
        else:
            raise TypeError("UFL objects not recognised.")

    form_data = tuple(_analyze_form(form, scalar_type) for form in forms)
    for data in form_data:
        elements += data.unique_sub_elements
        coordinate_elements += data.coordinate_elements

    # Loop through forms to extract interpolate operands
    new_coefficients = []
    for data in form_data:
        for interpolation in _interpolations(data):
            # Expose expression used for interpolation to generated code
            original_expression = interpolated_expression(interpolation)
            target_space = interpolation_target_space(interpolation)
            element = target_space.ufl_element()
            if isinstance(interpolation, ProxyCoefficient):
                # Map the expression onto the reference cell of the target
                # element, which is what its dual basis evaluates. An
                # interpolation that stayed in the integrands was mapped there
                # by `apply_interpolate_pullbacks`.
                domain = target_space.ufl_domain()
                mapped_expression = apply_push_forward(
                    original_expression, domain, element.pullback
                )
            else:
                mapped_expression = original_expression
            elements += ufl.algorithms.extract_elements(mapped_expression)
            processed_expression = _analyze_expression(
                mapped_expression,
                scalar_type,
                do_apply_function_pullbacks=isinstance(interpolation, ProxyCoefficient),
            )
            points = element.basix_element.points
            processed_expressions += [(processed_expression, points, original_expression)]

            # Append coefficents in the processed expression
            # to the form data reduced coefficients.
            new_coefficients.extend(
                [
                    coeff
                    for coeff in ufl.algorithms.extract_coefficients(processed_expression)
                    if coeff not in data.reduced_coefficients
                ]
            )
        data._reduced_coefficients.extend(new_coefficients)

        # Update form data for new set of reduced coefficients
        # Sort the reduced coefficients of the form.
        data._reduced_coefficients = sorted(data._reduced_coefficients, key=lambda x: x.count())
        # NOTE: Could have been simpler if we had extracted the elements from
        # reduced_coefficients rather than going through coefficient_elements.
        new_coeff_elements = tuple(coeff.ufl_element() for coeff in data.reduced_coefficients)
        data._coefficient_elements = new_coeff_elements
        # Enable all coefficients that are either in the original integral coefficients or
        # in the new coefficients generated from the interpolate expressions.
        for itg_data in data.integral_data:
            assert itg_data.integral_coefficients is not None
            itg_data.enabled_coefficients = [
                bool(coeff in itg_data.integral_coefficients) or bool(coeff in new_coefficients)
                for coeff in data.reduced_coefficients
            ]

        # Update original coefficient position in form
        data._original_coefficient_positions = [
            i
            for i, c in enumerate(data.original_form.coefficients())
            if c in data.reduced_coefficients
        ]

    for original_expression, points in expressions:
        elements += ufl.algorithms.extract_elements(original_expression)
        processed_expression = _analyze_expression(original_expression, scalar_type)
        processed_expressions += [(processed_expression, points, original_expression)]

    elements += ufl.algorithms.analysis.extract_sub_elements(elements)

    # Sort elements so sub-elements come before mixed elements
    unique_elements = ufl.algorithms.sort_elements(set(elements))
    unique_coordinate_element_list = sorted(set(coordinate_elements), key=lambda x: repr(x))

    for e in unique_elements:
        assert isinstance(e, basix.ufl._ElementBase)

    # Compute dict (map) from element to index
    element_numbers = {element: i for i, element in enumerate(unique_elements)}

    return UFLData(
        form_data=form_data,
        unique_elements=unique_elements,
        element_numbers=element_numbers,
        unique_coordinate_elements=unique_coordinate_element_list,
        expressions=processed_expressions,
    )


def _analyze_expression(
    expression: ufl.core.expr.Expr,
    scalar_type: npt.DTypeLike,
    do_apply_function_pullbacks: bool = True,
) -> ufl.core.expr.Expr:
    """Analyzes and preprocesses expressions.

    Args:
        expression: The expression to process.
        scalar_type: The scalar type to compile for.
        do_apply_function_pullbacks: Represent the form arguments and
            coefficients in reference value. This is already done for the
            operand of an interpolation that stayed in a form, which is
            collected after the form has been processed.
    """
    preserve_geometry_types = (ufl.classes.Jacobian,)
    expression = ufl.algorithms.apply_algebra_lowering.apply_algebra_lowering(expression)
    expression = ufl.algorithms.apply_derivatives.apply_derivatives(expression)
    if do_apply_function_pullbacks:
        expression = ufl.algorithms.apply_function_pullbacks.apply_function_pullbacks(expression)
    expression = ufl.algorithms.apply_geometry_lowering.apply_geometry_lowering(
        expression, preserve_geometry_types
    )
    expression = ufl.algorithms.apply_derivatives.apply_derivatives(expression)
    expression = ufl.algorithms.apply_geometry_lowering.apply_geometry_lowering(
        expression, preserve_geometry_types
    )
    expression = ufl.algorithms.apply_derivatives.apply_derivatives(expression)

    # Remove complex nodes if scalar type is real valued
    if not np.issubdtype(scalar_type, np.complexfloating):
        expression = ufl.algorithms.remove_complex_nodes.remove_complex_nodes(expression)

    return expression


def _analyze_form(form: ufl.Form, scalar_type: npt.DTypeLike) -> FormData:
    """Analyzes UFL form and attaches metadata.

    Args:
        form: forms
        scalar_type: Scalar type used for form. This is used to simplify
            real valued forms.

    Returns:
        Form data computed by UFL with metadata attached

    Note:
        The main workload of this function is extraction of
        unique/default metadata from options, integral metadata or
        inherited from UFL (in case of quadrature degree).
    """
    if form.empty():
        raise RuntimeError(f"Form ({form}) seems to be zero: cannot compile it.")
    if _has_custom_integrals(form):
        raise RuntimeError(f"Form ({form}) contains unsupported custom integrals.")

    # Check that coordinate element is based on basix.ufl._ElementBase
    for _integral in form._integrals:
        assert isinstance(_integral._ufl_domain._ufl_coordinate_element, basix.ufl._ElementBase)

    # Check for complex mode
    complex_mode = np.issubdtype(scalar_type, np.complexfloating)

    # Compute form metadata
    form_data: FormData = compute_form_data(
        form,
        do_apply_function_pullbacks=True,
        do_apply_integral_scaling=True,
        do_apply_geometry_lowering=True,
        preserve_geometry_types=(ufl.geometry.Jacobian,),  # type: ignore
        do_apply_restrictions=True,
        do_append_everywhere_integrals=False,  # do not add dx integrals to dx(i) in UFL
        complex_mode=complex_mode,
    )

    # Store original form, prior to replacement of interpolate,
    # to ensure we get the correct signature and coefficient numbering
    form_data._original_form = form

    # Determine unique quadrature degree and quadrature scheme
    # per each integral data
    for id, integral_data in enumerate(form_data.integral_data):
        # Iterate through groups of integral data. There is one integral
        # data for all integrals with same domain, itype, subdomain_id
        # (but possibly different metadata).
        #
        # Quadrature degree and quadrature scheme must be the same for
        # all integrals in this integral data group, i.e. must be the
        # same for for the same (domain, itype, subdomain_id)

        for i, integral in enumerate(integral_data.integrals):
            metadata = integral.metadata()

            # Vertex integrals do not support discontinuous integrands.
            if integral.integral_type() == "vertex":
                elements = ufl.algorithms.extract_elements(integral)
                if any(e.discontinuous for e in elements):
                    raise TypeError("Vertex integrals not supported for discontinuous elements.")

            # If form contains a quadrature element, use the custom
            # quadrature scheme
            custom_q = None
            for e in ufl.algorithms.extract_elements(integral):
                if e.has_custom_quadrature:
                    if custom_q is None:
                        custom_q = e.custom_quadrature()
                    else:
                        p, w = e.custom_quadrature()
                        assert np.allclose(p, custom_q[0])
                        assert np.allclose(w, custom_q[1])

            if custom_q is None:
                # Extract quadrature degree
                qd = -1
                if "quadrature_degree" in metadata.keys():
                    qd = metadata["quadrature_degree"]

                # Sending in a negative quadrature degree means that we want to be
                # able to customize it at a later stage.
                if qd < 0:
                    qd = int(np.max(integral.metadata()["estimated_polynomial_degree"]))
                # Extract quadrature rule
                qr = integral.metadata().get("quadrature_rule", "default")

                logger.info(f"Integral {i}, integral group {id}:")
                logger.info(f"--- quadrature rule: {qr}")
                logger.info(f"--- quadrature degree: {qd}")

                metadata.update({"quadrature_degree": qd, "quadrature_rule": qr})
            else:
                metadata.update(
                    {
                        "quadrature_points": custom_q[0],
                        "quadrature_weights": custom_q[1],
                        "quadrature_rule": "custom",
                    }
                )

            integral_data.integrals[i] = integral.reconstruct(metadata=metadata)

    return form_data


def _has_custom_integrals(
    o: ufl.integral.Integral | ufl.classes.Form | list | tuple,
) -> bool:
    """Check for custom integrals."""
    if isinstance(o, ufl.integral.Integral):
        return o.integral_type() in ufl.custom_integral_types
    elif isinstance(o, ufl.classes.Form):
        return any(_has_custom_integrals(itg) for itg in o.integrals())
    elif isinstance(o, list | tuple):
        return any(_has_custom_integrals(itg) for itg in o)
    else:
        raise NotImplementedError


class ProxyCoefficient(ufl.Coefficient):
    """Proxy coefficient to replace operands that require custom treatement."""

    _operator: ufl.core.expr.Expr

    def __init__(self, V: ufl.FunctionSpace, operator: ufl.core.expr.Expr):
        """Initialise."""
        self._operator = operator
        super().__init__(V)

    @property
    def operand(self) -> ufl.core.expr.Expr:
        """The operand that this proxy coefficient is replacing."""
        return self._operator.ufl_operands[0]

    @property
    def operator(self) -> ufl.core.expr.Expr:
        """The kind of the operand that this proxy coefficient is replacing."""
        return self._operator


def _primal_base_form_operator(
    o: ufl.core.base_form_operator.BaseFormOperator, argument: ufl.Argument
) -> ufl.core.base_form_operator.BaseFormOperator:
    """Rebuild a base form operator as a function of the space it is acting on.

    The operator taking part in an action does not know its own target space:
    its dual slot holds the form being acted on, so ``ufl_function_space`` is
    the space of the action's result. The space to rebuild it in is that of the
    argument the action contracts.

    Args:
        o: The base form operator to rebuild.
        argument: The argument that the action contracts.
    """
    (operand,) = o.ufl_operands
    return o._ufl_expr_reconstruct_(operand, v=argument.ufl_function_space())


def expand_base_form_operator_actions(form: ufl.form.BaseForm) -> ufl.Form:
    """Expand actions of base form operators into a plain form.

    Differentiating a form that holds an interpolation gives
    ``dF/dw = dF/dw|_N + sum_i Action(dF/dN_i, dN_i/dw)``, where ``N_i`` are the
    interpolations and ``dN_i/dw`` are two-forms. UFL folds an action into the
    operator's dual argument slot when it can, so the result is a
    {py:class}`ufl.form.FormSum` over forms, base form operators and actions.

    FFCx evaluates an interpolation cell-locally, so an action of it is not an
    operator that has to be assembled separately: substituting the interpolation
    for the argument that the left form acts on recovers an ordinary form.

    Args:
        form: The result of applying derivatives to a form.

    Return:
        An equivalent {py:class}`ufl.Form`.
    """
    if isinstance(form, ufl.Form):
        return form
    elif isinstance(form, ufl.form.FormSum):
        return sum(
            (
                weight * expand_base_form_operator_actions(component)
                for component, weight in zip(form.components(), form.weights())
            ),
            start=ufl.form.Form([]),
        )
    elif isinstance(form, ufl.classes.Action):
        left, right = form.ufl_operands
        if not isinstance(right, ufl.core.base_form_operator.BaseFormOperator):
            raise NotImplementedError(f"Cannot expand the action of {type(right).__name__}.")
        # `action` contracts the last argument of the left operand.
        left = expand_base_form_operator_actions(left)
        argument = left.arguments()[-1]
        return ufl.replace(left, {argument: _primal_base_form_operator(right, argument)})
    elif isinstance(form, ufl.core.base_form_operator.BaseFormOperator):
        # `Action(dF/dN, dN/dw)` folded into the dual argument slot, which is
        # only done for a single-argument left operand, so that argument is the
        # one standing in for the operator.
        dFdN = form.argument_slots()[0]
        dFdN = expand_base_form_operator_actions(dFdN)
        (argument,) = dFdN.arguments()
        return ufl.replace(dFdN, {argument: _primal_base_form_operator(form, argument)})
    else:
        raise NotImplementedError(f"Cannot expand a form of type {type(form).__name__}.")


class IntermediateCoefficientReplacer(DAGTraverser):
    """Replace operands requiring intermediate coefficients with intermediate objects."""

    def __init__(
        self,
        compress: bool | None = True,
        visited_cache: dict[tuple, ufl.core.expr.Expr] | None = None,
        result_cache: dict[ufl.core.expr.Expr, ufl.core.expr.Expr] | None = None,
    ) -> None:
        """Initialise.

        Args:
            compress: If True, ``result_cache`` will be used.
            visited_cache: cache of intermediate results;
                expr -> r = self.process(expr, ...).
            result_cache: cache of result objects for memory reuse, r -> r.

        """
        super().__init__(compress=compress, visited_cache=visited_cache, result_cache=result_cache)

    @singledispatchmethod
    def process(
        self,
        o: ufl.core.expr.Expr,
        reference_value: bool | None = False,
        reference_grad: int | None = 0,
        restricted: str | None = None,
    ) -> ufl.core.expr.Expr:
        """Replace averaged arguments with intermediate space.

        Args:
            o: `ufl.core.expr.Expr` to be processed.
            reference_value: Whether `ReferenceValue` has been applied or not.
            reference_grad: Number of `ReferenceGrad`s that have been applied.
            restricted: '+', '-', or None.
        """
        return super().process(o)

    @process.register(ufl.Interpolate)
    def _(
        self,
        o: ufl.Interpolate,
        reference_value: bool | None = False,
        reference_grad: int | None = 0,
        restricted: str | None = None,
    ) -> ufl.core.expr.Expr:
        """Handle Interpolate."""
        (operand,) = o.ufl_operands
        if ufl.algorithms.extract_arguments(operand):
            # An interpolation holding an argument stays in the integrand: it is
            # a modified terminal, see `ffcx.ir.analysis.modified_terminals`.
            return self.reuse_if_untouched(
                o,
                reference_value=reference_value,
                reference_grad=reference_grad,
                restricted=restricted,
            )
        return ProxyCoefficient(interpolation_target_space(o), o)

    @process.register(ufl.core.expr.Expr)
    def _(
        self,
        o: ufl.Argument,
        reference_value: bool | None = False,
        reference_grad: int | None = 0,
        restricted: str | None = None,
    ) -> ufl.core.expr.Expr:
        """Handle anything else in UFL."""
        return self.reuse_if_untouched(
            o,
            reference_value=reference_value,
            reference_grad=reference_grad,
            restricted=restricted,
        )


def replace_ufl_operands(form: ufl.Form) -> ufl.Form:
    """Parse UFL form and and add placeholder coefficients.

    Args:
        form: UFL form

    Return:
        The modified form with operands replaced.
    """
    rule = IntermediateCoefficientReplacer()
    return ufl.algorithms.map_integrands.map_integrands(rule, form)  # type: ignore


def compute_form_data(
    form: ufl.Form,
    do_apply_function_pullbacks: bool = False,
    do_apply_integral_scaling: bool = False,
    do_apply_geometry_lowering: bool = False,
    preserve_geometry_types: tuple[ufl.geometry.GeometricQuantity, ...] = (),
    do_apply_default_restrictions: bool = True,
    do_apply_restrictions: bool = True,
    do_estimate_degrees: bool = True,
    do_append_everywhere_integrals: bool = True,
    do_replace_functions: bool = False,
    coefficients_to_split: tuple[ufl.Coefficient, ...] | None = None,
    complex_mode: bool = False,
    do_remove_component_tensors: bool = False,
) -> FormData:
    """Compute form data.

    Args:
        form: The form to compute form data for.
        do_apply_function_pullbacks: Apply pull-back to reference cell
            for coefficients, including Piola and symmetry transforms
            if required.
        do_apply_integral_scaling: Apply scaling of moving the integral
            from physical to reference frame.
        do_apply_geometry_lowering: Lower the representation of geometrical
            quantities to a smaller subset of quantities
        preserve_geometry_types: Set of quantities not to lower, and keep
            at its present stage for the form-compiler.
        do_apply_default_restrictions: Apply default restrictions, defined in
            {py:mod}`ufl.algorithms.apply_restrictions` to integrals if no
            restriction has been set.
        do_apply_restrictions: Apply restrictions towards terminal nodes.
        do_replace_functions: Replace functions with with its cannonically numbered
            function or thos provided in coefficients_to_split.
        coefficients_to_split: Sequence of coefficients to split over a MeshSequence.
        do_estimate_degrees: Estimate polynomial degree of integrands.
        do_append_everywhere_integrals: If True append every `dx` integral to each `dx(i)`
            integral defined in the form.
        do_remove_component_tensors: Remove component-tensor if true.
        complex_mode: If false remove complex nodes from the form.
    """
    # --- Store untouched form for reference.
    # The user of FormData may get original arguments,
    # original coefficients, and form signature from this object.
    # But be aware that the set of original coefficients are not
    # the same as the ones used in the final UFC form.
    # See 'reduced_coefficients' below.

    # --- Expand derivatives of base form operators.
    # An interpolation is linear in the expression it interpolates, so a
    # derivative can be pushed into it: d/du I(g(u))[v] = I(dg/du[v]). That
    # keeps the result an ordinary form and composes to higher derivatives,
    # which the Jacobian of an interpolated expression needs. `preprocess_form`
    # below repeats both steps, which is idempotent.
    form = apply_derivatives(apply_algebra_lowering(form), differentiate_through_operators=True)

    # --- Expand any remaining action of a base form operator into a plain form.
    form = expand_base_form_operator_actions(form)

    original_form = form

    # --- Pass form integrands through some symbolic manipulation
    # Interpolations of an expression become proxy coefficients. These are
    # internal to the generated kernel, so they must *not* appear in
    # `original_form`, whose coefficients are the ones the caller passes in.
    form = replace_ufl_operands(form)

    # Evaluate the interpolations that are left, which are the ones holding an
    # argument, on the reference cell of their target element. The ones that
    # became proxy coefficients are gone by now: their operand is compiled as an
    # expression of its own, and mapped there in `analyze_ufl_objects`.
    if do_apply_function_pullbacks:
        form = apply_interpolate_pullbacks(form)
        for integral in form.integrals():
            for node in ufl.corealg.traversal.unique_pre_traversal(integral.integrand()):
                if isinstance(node, ufl.Interpolate):
                    check_interpolation(node)

    form = preprocess_form(form, complex_mode)

    # --- Group form integrals
    # TODO: Refactor this, it's rather opaque what this does
    # TODO: Is self.original_form.ufl_domains() right here?
    #       It will matter when we start including 'num_domains' in ufc form.
    form = group_form_integrals(
        form,
        original_form.ufl_domains(),
        do_append_everywhere_integrals=do_append_everywhere_integrals,
    )

    # Estimate polynomial degree of integrands now, before applying
    # any pullbacks and geometric lowering.  Otherwise quad degrees
    # blow up horrifically.
    if do_estimate_degrees:
        form = attach_estimated_degrees(form)

    if do_apply_function_pullbacks:
        # Rewrite coefficients and arguments in terms of their
        # reference cell values with Piola transforms and symmetry
        # transforms injected where needed.
        # Decision: Not supporting grad(dolfin.Expression) without a
        #           Domain.  Current dolfin works if Expression has a
        #           cell but this should be changed to a mesh.
        form = apply_function_pullbacks(form)

    # Scale integrals to reference cell frames
    if do_apply_integral_scaling:
        form = apply_integral_scaling(form)

    # Lower abstractions for geometric quantities into a smaller set
    # of quantities, allowing the form compiler to deal with a smaller
    # set of types and treating geometric quantities like any other
    # expressions w.r.t. loop-invariant code motion etc.
    if do_apply_geometry_lowering:
        form = apply_geometry_lowering(form, preserve_geometry_types)

    # Apply differentiation again, because the algorithms above can
    # generate new derivatives or rewrite expressions inside
    # derivatives
    if do_apply_function_pullbacks or do_apply_geometry_lowering:
        form = apply_derivatives(form)

        # Neverending story: apply_derivatives introduces new Jinvs,
        # which needs more geometry lowering
        if do_apply_geometry_lowering:
            form = apply_geometry_lowering(form, preserve_geometry_types)
            # Lower derivatives that may have appeared
            form = apply_derivatives(form)

    form = apply_coordinate_derivatives(form)

    # If in real mode, remove any complex nodes introduced during form processing.
    if not complex_mode:
        form = remove_complex_nodes(form)

    # Remove component tensors
    if do_remove_component_tensors:
        form = remove_component_tensors(form)
    integral_data = build_integral_data(form.integrals())
    return FormData(
        original_form,
        integral_data,
        do_apply_default_restrictions=do_apply_default_restrictions,
        do_apply_restrictions=do_apply_restrictions,
        do_replace_functions=do_replace_functions,
        coefficients_to_split=coefficients_to_split,
        complex_mode=complex_mode,
    )
