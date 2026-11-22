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
from ufl.algorithms.apply_derivatives import apply_derivatives
from ufl.algorithms.apply_function_pullbacks import apply_interpolate_pullbacks
from ufl.algorithms.apply_geometry_lowering import apply_geometry_lowering
from ufl.algorithms.remove_component_tensors import remove_component_tensors
from ufl.corealg.dag_traverser import DAGTraverser
from ufl.corealg.traversal import unique_pre_traversal
from ufl.utils.sorting import sorted_by_count

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
    #: The expression kernel that evaluates each interpolation at the points
    #: of its target element, as an index into `expressions`
    interpolation_kernels: dict[ProxyCoefficient | ufl.Interpolate, int]


def kernel_interpolations(
    expressions: typing.Iterable[ufl.core.expr.Expr],
    coefficients: typing.Iterable[ufl.Coefficient],
) -> list[ProxyCoefficient | ufl.Interpolate]:
    """Collect what needs an expression kernel evaluated at interpolation points.

    A proxy coefficient always needs one: its degrees of freedom are built per
    cell by evaluating the interpolated expression. An interpolation that
    survived into the expressions, which is one holding an argument, needs one
    too, unless its expression is the argument itself and the table is known at
    compile time.

    Args:
        expressions: Processed expressions, such as the integrands of a form.
        coefficients: The coefficients of the expressions.
    """
    collected: list[ProxyCoefficient | ufl.Interpolate] = [
        coefficient for coefficient in coefficients if isinstance(coefficient, ProxyCoefficient)
    ]
    seen = set()
    for expression in expressions:
        for node in unique_pre_traversal(expression):
            if (
                isinstance(node, ufl.Interpolate)
                and node not in seen
                and interpolation_has_runtime_table(node)
            ):
                seen.add(node)
                collected.append(node)
    return collected


def required_coefficients(expression: ufl.core.expr.Expr) -> list[ufl.Coefficient]:
    """Get the coefficients that evaluating an expression reads, sorted by count.

    A proxy coefficient is built from the coefficients of the expression it
    interpolates, so it stands for those rather than for itself.

    Args:
        expression: An expression.
    """
    coefficients: set[ufl.Coefficient] = set()
    for coefficient in ufl.algorithms.extract_coefficients(expression):
        if isinstance(coefficient, ProxyCoefficient):
            coefficients.update(required_coefficients(coefficient.operand))
        else:
            coefficients.add(coefficient)
    return typing.cast(list[ufl.Coefficient], sorted_by_count(coefficients))


def required_constants(expression: ufl.core.expr.Expr) -> list[ufl.Constant]:
    """Get the constants that evaluating an expression reads, sorted by count.

    Args:
        expression: An expression.
    """
    constants = set(ufl.algorithms.analysis.extract_constants(expression))
    for coefficient in ufl.algorithms.extract_coefficients(expression):
        if isinstance(coefficient, ProxyCoefficient):
            constants.update(required_constants(coefficient.operand))
    return typing.cast(list[ufl.Constant], sorted_by_count(constants))


def interpolation_coefficients(
    interpolation: ProxyCoefficient | ufl.Interpolate,
) -> list[ufl.Coefficient]:
    """Get the coefficients that the expression kernel of an interpolation reads.

    They are passed to the kernel in this order. It is taken from the expression
    as given rather than as compiled, which may have lost some, so that two
    interpolations that differ only in their coefficients share one kernel.

    Args:
        interpolation: A `ProxyCoefficient` or a `ufl.Interpolate`.
    """
    return required_coefficients(interpolated_expression(interpolation))


def interpolation_constants(
    interpolation: ProxyCoefficient | ufl.Interpolate,
) -> list[ufl.Constant]:
    """Get the constants that the expression kernel of an interpolation reads.

    They are passed to the kernel in this order, see `interpolation_coefficients`.

    Args:
        interpolation: A `ProxyCoefficient` or a `ufl.Interpolate`.
    """
    return required_constants(interpolated_expression(interpolation))


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
    # reference value, through its own pull back. Both forms are met: an
    # expression is checked between the two passes, and a form after them.
    # Cancelling Jacobian products in a form then removes the two pull backs
    # entirely, which `check_interpolation` requires to be the same.
    mapped = (
        target.pullback.apply_inverse(argument, domain),
        ufl.classes.ReferenceValue(argument),
        target.pullback.apply_inverse(
            source.pullback.apply(ufl.classes.ReferenceValue(argument), domain), domain
        ),
    )
    return argument if any(_same_expression(operand, m) for m in mapped) else None


def _same_expression(a: ufl.core.expr.Expr, b: ufl.core.expr.Expr) -> bool:
    """Compare two expressions up to lowering and index numbering.

    The operand of an interpolation is compared against one built here, which
    has not been through the passes the form has: the geometry of a pull back is
    lowered in the form, cancelling Jacobian products removes its component
    tensors, and each pull back numbers its free indices from a global counter.

    Args:
        a: An expression.
        b: An expression to compare it against.
    """

    def normalise(expression: ufl.core.expr.Expr) -> ufl.core.expr.Expr:
        preserve = (ufl.classes.Jacobian,)
        for _ in range(2):
            expression = apply_geometry_lowering(expression, preserve)
            expression = apply_derivatives(expression)
        expression = remove_component_tensors(expression)
        return typing.cast(
            ufl.core.expr.Expr, ufl.algorithms.renumbering.renumber_indices(expression)
        )

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
        return typing.cast(ufl.FunctionSpace, interpolation.ufl_function_space())
    # `ufl_function_space()` is the dual of the source space for an interpolation
    # of a test function, so build the target space from the target element.
    domain = ufl.domain.extract_unique_domain(interpolation)
    return ufl.FunctionSpace(domain, interpolation.ufl_element())


def analyze_ufl_objects(
    ufl_objects: list[
        ufl.form.Form
        | basix.ufl._ElementBase
        | ufl.Mesh
        | tuple[ufl.core.expr.Expr, npt.NDArray[np.floating]]
    ],
    scalar_type: npt.DTypeLike,
    # do_cancel_jacobian_products: bool = True,
    do_cancel_jacobian_products: bool = False,
) -> UFLData:
    """Analyze ufl object(s).

    Args:
        ufl_objects: UFL objects
        scalar_type: Scalar type that should be used for the analysis
        do_cancel_jacobian_products: Whether to cancel jacobian products.

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

    form_data = tuple(
        _analyze_form(form, scalar_type, do_cancel_jacobian_products) for form in forms
    )
    for data in form_data:
        elements += data.unique_sub_elements
        coordinate_elements += data.coordinate_elements

    # Each interpolation that needs an expression kernel has it compiled as an
    # expression of its own. A kernel is added before the expressions that call
    # it, including the kernels of interpolations nested in its own expression,
    # so that it is generated first.
    interpolation_kernels: dict[ProxyCoefficient | ufl.Interpolate, int] = {}

    def add_interpolation_kernels(
        interpolations: typing.Iterable[ProxyCoefficient | ufl.Interpolate],
    ) -> None:
        for interpolation in interpolations:
            if interpolation in interpolation_kernels:
                continue
            processed, points, mapped = _analyze_interpolation(interpolation, scalar_type)
            elements.extend(ufl.algorithms.extract_elements(mapped))
            elements.extend(ufl.algorithms.extract_elements(processed))
            add_interpolation_kernels(
                kernel_interpolations([processed], ufl.algorithms.extract_coefficients(processed))
            )
            interpolation_kernels[interpolation] = len(processed_expressions)
            processed_expressions.append((processed, points, mapped))

    for data in form_data:
        # The coefficients that the interpolations in each integral read must be
        # passed to it, even if the integrand does not hold them itself.
        interpolated_coefficients: list[set[ufl.Coefficient]] = []
        for itg_data in data.integral_data:
            assert itg_data.integral_coefficients is not None
            interpolations = kernel_interpolations(
                [integral.integrand() for integral in itg_data.integrals],
                [c for c in data.reduced_coefficients if c in itg_data.integral_coefficients],
            )
            add_interpolation_kernels(interpolations)
            interpolated_coefficients.append(
                {c for i in interpolations for c in interpolation_coefficients(i)}
            )

        # Update form data for new set of reduced coefficients
        data._reduced_coefficients = sorted_by_count(
            set(data.reduced_coefficients).union(*interpolated_coefficients)
        )
        # NOTE: Could have been simpler if we had extracted the elements from
        # reduced_coefficients rather than going through coefficient_elements.
        data._coefficient_elements = tuple(c.ufl_element() for c in data.reduced_coefficients)
        for itg_data, interpolated in zip(data.integral_data, interpolated_coefficients):
            assert itg_data.integral_coefficients is not None
            itg_data.enabled_coefficients = [
                coeff in itg_data.integral_coefficients or coeff in interpolated
                for coeff in data.reduced_coefficients
            ]

        # Update original coefficient position in form
        data._original_coefficient_positions = [
            i
            for i, c in enumerate(data.original_form.coefficients())
            if c in data.reduced_coefficients
        ]

    for original_expression, points in expressions:
        processed_expression = _analyze_expression(original_expression, scalar_type)
        # Taken from the processed expression: in the original, UFL counts the
        # coargument of an interpolation that is not inside another expression
        # as a free argument, which clashes with those of nested ones.
        elements += ufl.algorithms.extract_elements(processed_expression)
        add_interpolation_kernels(
            kernel_interpolations(
                [processed_expression], ufl.algorithms.extract_coefficients(processed_expression)
            )
        )
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
        interpolation_kernels=interpolation_kernels,
    )


def _analyze_interpolation(
    interpolation: ProxyCoefficient | ufl.Interpolate, scalar_type: npt.DTypeLike
) -> tuple[ufl.core.expr.Expr, npt.NDArray[np.floating], ufl.core.expr.Expr]:
    """Analyze the expression that an interpolation evaluates at interpolation points.

    Args:
        interpolation: A `ProxyCoefficient` or a `ufl.Interpolate`.
        scalar_type: The scalar type to compile for.

    Returns:
        The processed expression, the interpolation points of the target
        element, and the expression on the reference cell of the target element.
    """
    expression = interpolated_expression(interpolation)
    target_space = interpolation_target_space(interpolation)
    element = target_space.ufl_element()
    if isinstance(interpolation, ProxyCoefficient):
        # Map the expression onto the reference cell of the target element,
        # which is what its dual basis evaluates. An interpolation that stayed
        # in the integrands was mapped there by `apply_interpolate_pullbacks`.
        expression = apply_push_forward(expression, target_space.ufl_domain(), element.pullback)
    processed = _analyze_expression(
        expression,
        scalar_type,
        do_apply_function_pullbacks=isinstance(interpolation, ProxyCoefficient),
    )
    return processed, element.basix_element.points, expression


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
    lowered = ufl.algorithms.apply_algebra_lowering.apply_algebra_lowering(expression)
    assert isinstance(lowered, ufl.core.expr.Expr)
    # An interpolation is linear in the expression it interpolates, so a
    # derivative is pushed into it, see `compute_form_data`.
    expression = ufl.algorithms.apply_derivatives.apply_derivatives(lowered)
    # Interpolations of an expression become proxy coefficients, which are
    # evaluated by an expression kernel of their own.
    expression = typing.cast(ufl.core.expr.Expr, IntermediateCoefficientReplacer()(expression))
    if do_apply_function_pullbacks:
        # Evaluate the interpolations that are left, which are the ones holding
        # an argument, on the reference cell of their target element.
        expression = apply_interpolate_pullbacks(expression)
        for node in unique_pre_traversal(expression):
            if isinstance(node, ufl.Interpolate):
                check_interpolation(node)
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


def _analyze_form(
    form: ufl.Form, scalar_type: npt.DTypeLike, do_cancel_jacobian_products: bool
) -> FormData:
    """Analyzes UFL form and attaches metadata.

    Args:
        form: forms
        scalar_type: Scalar type used for form. This is used to simplify
            real valued forms.
        do_cancel_jacobian_products: Whether to cancel jacobian products.

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
        preserve_geometry_types=(ufl.classes.Jacobian,),
        do_apply_restrictions=True,
        do_append_everywhere_integrals=False,  # do not add dx integrals to dx(i) in UFL
        do_cancel_jacobian_products=do_cancel_jacobian_products,
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

            # Extract from the integrand, an expression, so that the
            # Coargument of an interpolation in it counts as bound, as it is
            # in the form. Passed the integral, UFL takes it for a free
            # argument, which clashes with the form's own argument of that
            # number once a derivative has been pushed into the interpolation.
            elements = ufl.algorithms.extract_elements(integral.integrand())

            # Vertex integrals do not support discontinuous integrands.
            if integral.integral_type() == "vertex":
                if any(e.discontinuous for e in elements):
                    raise TypeError("Vertex integrals not supported for discontinuous elements.")

            # If form contains a quadrature element, use the custom
            # quadrature scheme
            custom_q = None
            for e in elements:
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

                # An upper bound on the estimated degree. It caps the
                # degree that UFL estimates, and is contradictory if the
                # degree has also been set explicitly.
                qd_max = metadata.get("max_quadrature_degree")
                if qd_max is not None:
                    if qd >= 0:
                        raise ValueError(
                            "Cannot set both 'quadrature_degree' and "
                            "'max_quadrature_degree' in integral metadata. "
                            "'max_quadrature_degree' bounds the estimated degree, "
                            "so it has no effect when the degree is given explicitly."
                        )
                    if int(qd_max) < 0:
                        raise ValueError(
                            f"'max_quadrature_degree' must be non-negative, got {qd_max}."
                        )

                # Sending in a negative quadrature degree means that we want to be
                # able to customize it at a later stage.
                if qd < 0:
                    qd = int(np.max(integral.metadata()["estimated_polynomial_degree"]))
                    if qd_max is not None:
                        qd = min(qd, int(qd_max))
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


class IntermediateCoefficientReplacer(DAGTraverser):
    """Replace operands requiring intermediate coefficients with intermediate objects."""

    def __init__(
        self,
        compress: bool | None = True,
        visited_cache: dict[tuple, ufl.core.expr.Expr | ufl.form.BaseForm] | None = None,
        result_cache: dict[
            ufl.core.expr.Expr | ufl.form.BaseForm, ufl.core.expr.Expr | ufl.form.BaseForm
        ]
        | None = None,
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
            reference_grad: Number of `ReferenceGrad` operators that have been applied.
            restricted: '+', '-', or None.
        """
        return typing.cast(ufl.core.expr.Expr, super().process(o))

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


def compute_form_data(form: ufl.Form, **kwargs: typing.Any) -> FormData:
    """Compute form data, with interpolations of expressions as proxy coefficients.

    Wraps {py:func}`ufl.algorithms.compute_form_data`, which it passes
    ``kwargs`` on to, after replacing each interpolation that holds no argument
    by a `ProxyCoefficient`. The remaining interpolations, which hold an
    argument, are evaluated on the reference cell of their target element by
    UFL, and checked here.

    Note:
        The `FormData.original_form` returned holds the proxy coefficients.

    Args:
        form: The form to compute form data for.
        kwargs: Options of {py:func}`ufl.algorithms.compute_form_data`.
    """
    # An interpolation is linear in the expression it interpolates, so UFL
    # pushes a derivative into it: d/du I(g(u))[v] = I(dg/du[v]). This has to
    # happen before the interpolation is replaced, since a proxy coefficient is
    # opaque to differentiation. `compute_form_data` repeats both steps, which
    # is idempotent.
    form = apply_derivatives(apply_algebra_lowering(form))
    form = replace_ufl_operands(form)

    form_data = ufl.algorithms.compute_form_data(form, **kwargs)

    if kwargs.get("do_apply_function_pullbacks", False):
        for itg_data in form_data.integral_data:
            for integral in itg_data.integrals:
                for node in unique_pre_traversal(integral.integrand()):
                    if isinstance(node, ufl.Interpolate):
                        check_interpolation(node)
    return form_data
