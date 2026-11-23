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
from ufl.algorithms.apply_derivatives import apply_derivatives
from ufl.algorithms.apply_function_pullbacks import apply_interpolate_pullbacks
from ufl.algorithms.apply_geometry_lowering import apply_geometry_lowering
from ufl.algorithms.remove_component_tensors import remove_component_tensors
from ufl.corealg.dag_traverser import DAGTraverser
from ufl.corealg.traversal import unique_pre_traversal
from ufl.utils.sorting import sorted_by_count

logger = logging.getLogger("ffcx")

#: The pull backs of the spaces that a proxy coefficient can be interpolated into
_PROXY_PULLBACKS = (
    ufl.pullback.IdentityPullback,
    ufl.pullback.L2Piola,
    ufl.pullback.CovariantPiola,
    ufl.pullback.ContravariantPiola,
    ufl.pullback.DoubleCovariantPiola,
    ufl.pullback.DoubleContravariantPiola,
    ufl.pullback.CovariantContravariantPiola,
)


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
) -> list[ProxyCoefficient | ufl.Interpolate]:
    """Collect what needs an expression kernel evaluated at interpolation points.

    A proxy coefficient always needs one: its degrees of freedom are built per
    cell by evaluating the interpolated expression. An interpolation that
    survived into the expressions, which is one holding an argument, needs one
    too, unless its expression is the argument itself and the table is known at
    compile time.

    Args:
        expressions: Processed expressions, such as the integrands of a form.

    Returns:
        The proxy coefficients, ordered by count, followed by the interpolations
        whose table is built for each cell, in the order they are met.
    """
    expressions = list(expressions)
    proxies = {
        coefficient
        for expression in expressions
        for coefficient in ufl.algorithms.extract_coefficients(expression)
        if isinstance(coefficient, ProxyCoefficient)
    }
    collected: list[ProxyCoefficient | ufl.Interpolate] = list(sorted_by_count(proxies))
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


def integral_interpolations(
    itg_data: ufl.algorithms.domain_analysis.IntegralData,
) -> list[ProxyCoefficient | ufl.Interpolate]:
    """Collect what needs an expression kernel evaluated in an integral.

    See `kernel_interpolations`. The proxy coefficients come first, ordered by
    count, which is also the order they are numbered in within the kernel.

    Args:
        itg_data: The integral data of a form analysed by `_analyze_form`,
            whose integrands hold the proxy coefficients.
    """
    return kernel_interpolations(integral.integrand() for integral in itg_data.integrals)


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
    # `apply_function_pullbacks` then represents the argument itself in
    # reference value, through its own pull back. Interpolations are only
    # classified after both passes. Cancelling Jacobian products in a form
    # then removes the two pull backs entirely, which `check_interpolation`
    # requires to be the same.
    mapped = (
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


def check_proxy_targets(o: ufl.Form | ufl.core.expr.Expr) -> None:
    """Check that each interpolation of an expression lands in a supported space.

    Such an interpolation becomes a `ProxyCoefficient`, whose dofs are built by
    `ffcx.codegeneration.interpolation` from a single basix element, which a
    mixed, symmetric or quadrature space lacks. This is checked before UFL maps
    the interpolations onto their reference cells, which fails for some of them.

    Args:
        o: A form or an expression, as given.
    """
    expressions = [itg.integrand() for itg in o.integrals()] if isinstance(o, ufl.Form) else [o]
    for expression in expressions:
        for node in unique_pre_traversal(expression):
            if not isinstance(node, ufl.Interpolate) or interpolation_arguments(node):
                continue
            pullback = node.ufl_element().pullback
            if not isinstance(pullback, _PROXY_PULLBACKS):
                raise NotImplementedError(
                    f"Interpolating an expression into a {pullback} space is not supported."
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
    do_cancel_jacobian_products: bool = True,
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
            # A proxy coefficient is not a coefficient of the form, so its
            # element is not among those of the form data.
            target = interpolation.ufl_element()
            assert isinstance(target, basix.ufl._ElementBase)
            elements.append(target)
            elements.extend(ufl.algorithms.extract_elements(mapped))
            elements.extend(ufl.algorithms.extract_elements(processed))
            add_interpolation_kernels(kernel_interpolations([processed]))
            interpolation_kernels[interpolation] = len(processed_expressions)
            processed_expressions.append((processed, points, mapped))

    # The coefficients that the interpolations read are coefficients of the
    # form, which UFL collects from inside the interpolations it processed.
    for data in form_data:
        for itg_data in data.integral_data:
            add_interpolation_kernels(integral_interpolations(itg_data))

    for original_expression, points in expressions:
        processed_expression = _analyze_expression(original_expression, scalar_type)
        # Taken from the processed expression: in the original, UFL counts the
        # coargument of an interpolation that is not inside another expression
        # as a free argument, which clashes with those of nested ones.
        elements += ufl.algorithms.extract_elements(processed_expression)
        add_interpolation_kernels(kernel_interpolations([processed_expression]))
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
    # `apply_interpolate_pullbacks` has mapped the expression onto the reference
    # cell of the target element, which is what its dual basis evaluates, and
    # `apply_function_pullbacks` has represented what it holds in reference value.
    expression = interpolated_expression(interpolation)
    processed = _analyze_expression(expression, scalar_type, do_apply_function_pullbacks=False)
    element = interpolation.ufl_element()
    assert isinstance(element, basix.ufl._ElementBase)
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
        do_apply_function_pullbacks: Evaluate the interpolations on the
            reference cells of their target elements, and represent the form
            arguments and coefficients in reference value. This is already done
            for the expression an interpolation interpolates, which is taken
            from a processed form or expression.
    """
    preserve_geometry_types = (ufl.classes.Jacobian,)
    lowered = ufl.algorithms.apply_algebra_lowering.apply_algebra_lowering(expression)
    assert isinstance(lowered, ufl.core.expr.Expr)
    # An interpolation is linear in the expression it interpolates, so UFL
    # pushes a derivative into it.
    expression = ufl.algorithms.apply_derivatives.apply_derivatives(lowered)
    if do_apply_function_pullbacks:
        check_proxy_targets(expression)
        # The same passes, in the same order, as `ufl.algorithms.compute_form_data`
        expression = apply_interpolate_pullbacks(expression)
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

    return _replace_interpolations(expression)


def _replace_interpolations(
    expression: ufl.core.expr.Expr, replacer: IntermediateCoefficientReplacer | None = None
) -> ufl.core.expr.Expr:
    """Replace the interpolations of an expression by proxy coefficients.

    Each one is evaluated by an expression kernel of its own. The
    interpolations that are left hold an argument, and are checked to be
    supported.

    Args:
        expression: A processed expression, or integrand of a processed form.
        replacer: The replacer to use. One that is shared between expressions
            replaces an interpolation they share by the same proxy coefficient.
    """
    replacer = replacer or IntermediateCoefficientReplacer()
    expression = typing.cast(ufl.core.expr.Expr, replacer(expression))
    for node in unique_pre_traversal(expression):
        if isinstance(node, ufl.Interpolate):
            check_interpolation(node)
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

    # Interpolations of an expression are left in the form for UFL to process,
    # so that it accounts for the coefficients they read, and replaced by proxy
    # coefficients below. Fail early for the ones that cannot be.
    check_proxy_targets(form)

    # Compute form metadata
    form_data: FormData = ufl.algorithms.compute_form_data(
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

    # Shared between the integrals, so that an interpolation they share is
    # replaced by the same proxy coefficient
    replacer = IntermediateCoefficientReplacer()

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

            # Interpolations of an expression become proxy coefficients. UFL
            # has mapped each onto the reference cell of its target element,
            # and represented what it reads in reference value.
            integrand = _replace_interpolations(integral.integrand(), replacer)

            # Extract from the integrand, an expression, so that the
            # Coargument of an interpolation in it counts as bound, as it is
            # in the form. Passed the integral, UFL takes it for a free
            # argument, which clashes with the form's own argument of that
            # number once a derivative has been pushed into the interpolation.
            elements = ufl.algorithms.extract_elements(integrand)

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

            integral_data.integrals[i] = integral.reconstruct(
                integrand=integrand, metadata=metadata
            )

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
    """A coefficient that stands for an interpolation of an expression.

    Its degrees of freedom are built in the kernel, by evaluating the expression
    at the interpolation points of the target element, see
    `ffcx.codegeneration.interpolation`.
    """

    _operator: ufl.Interpolate

    def __init__(self, V: ufl.FunctionSpace, operator: ufl.Interpolate):
        """Initialise.

        Args:
            V: The space interpolated into.
            operator: The interpolation, of an expression already mapped onto
                the reference cell of the target element.
        """
        self._operator = operator
        super().__init__(V)

    @property
    def operand(self) -> ufl.core.expr.Expr:
        """The expression interpolated, on the reference cell of the target element."""
        return self._operator.ufl_operands[0]

    @property
    def operator(self) -> ufl.Interpolate:
        """The interpolation that this proxy coefficient stands for."""
        return self._operator


class IntermediateCoefficientReplacer(DAGTraverser):
    """Replace interpolations of an expression by proxy coefficients.

    An interpolation holding an argument is left in place: it is a modified
    terminal, see `ffcx.ir.analysis.modified_terminals`. Interpolations nested in
    the expression of a replaced one are left too, and replaced when that
    expression is compiled into a kernel of its own.
    """

    @singledispatchmethod
    def process(self, o: ufl.core.expr.Expr) -> ufl.core.expr.Expr:
        """Replace the interpolations of an expression.

        Args:
            o: The expression.
        """
        return typing.cast(ufl.core.expr.Expr, super().process(o))

    @process.register(ufl.Interpolate)
    def _(self, o: ufl.Interpolate) -> ufl.core.expr.Expr:
        """Handle an interpolation."""
        if interpolation_arguments(o):
            return self.reuse_if_untouched(o)
        return ProxyCoefficient(interpolation_target_space(o), o)

    @process.register(ufl.core.expr.Expr)
    def _(self, o: ufl.core.expr.Expr) -> ufl.core.expr.Expr:
        """Handle anything else."""
        return self.reuse_if_untouched(o)
