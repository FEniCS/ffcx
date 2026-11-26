# Copyright (C) 2007-2026 Anders Logg, Martin Alnaes, Kristian B. Oelgaard,
#                         Michal Habera, Jørgen S. Dokken, Paul T. Kühner, and others
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

import functools
import itertools
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
from ufl.algorithms.apply_geometry_lowering import apply_geometry_lowering
from ufl.algorithms.renumbering import renumber_indices
from ufl.corealg.dag_traverser import DAGTraverser
from ufl.corealg.traversal import unique_pre_traversal

logger = logging.getLogger("ffcx")


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


def interpolated_expression(interpolation: ufl.Interpolate) -> ufl.core.expr.Expr:
    """Get the expression that an interpolation interpolates.

    Args:
        interpolation: A `ufl.Interpolate`.
    """
    (operand,) = interpolation.ufl_operands
    return operand


def interpolation_arguments(interpolation: ufl.Interpolate) -> tuple[ufl.Argument, ...]:
    """Get the arguments that an interpolated expression is linear in, by number.

    Args:
        interpolation: A `ufl.Interpolate`.
    """
    arguments = ufl.algorithms.extract_arguments(interpolated_expression(interpolation))
    return tuple(sorted(arguments, key=lambda argument: argument.number()))


def interpolated_argument(interpolation: ufl.Interpolate) -> ufl.Argument | None:
    """Get the argument that an analysed interpolation maps directly, if it maps one.

    The analysis rewrites the operand of each interpolation of an argument
    itself into the argument's reference value, see `_classify_interpolation`,
    so after it this is an exact comparison.

    Args:
        interpolation: A `ufl.Interpolate`, from an analysed form or expression.

    Returns:
        The argument, or None if the operand is not its reference value.
    """
    arguments = interpolation_arguments(interpolation)
    if len(arguments) != 1:
        return None
    (argument,) = arguments
    if interpolated_expression(interpolation) != ufl.classes.ReferenceValue(argument):
        return None
    return argument


class _Processing(typing.NamedTuple):
    """How an expression holding interpolations has been processed.

    An interpolation in it is classified against an interpolation of the
    argument itself processed the same way, see `_processed_argument`.
    """

    #: As the integrand of a form, by `_compute_form_data`, rather than as an
    #: expression, by `_process_expression`
    form: bool
    #: For a complex scalar type
    complex_mode: bool
    #: With Jacobian products cancelled, which only applies to a form
    do_cancel_jacobian_products: bool = False


@functools.lru_cache(maxsize=128)
def _processed_argument(
    space: ufl.FunctionSpace, part: int | None, target: ufl.FunctionSpace, processing: _Processing
) -> tuple[ufl.Argument, ufl.core.expr.Expr]:
    """Get what interpolating an argument itself is processed into.

    The interpolation of an argument into the target space is processed the
    way that `processing` describes, so that it carries what UFL makes of it:
    the pull backs mapping it onto the reference cell of the target element,
    their lowered geometry and, if enabled, cancelled Jacobian products.

    Args:
        space: The space of the argument.
        part: The part of the argument.
        target: The space interpolated into.
        processing: How to process the interpolation.

    Returns:
        The argument it is built with, and the operand of the processed
        interpolation, with its indices renumbered.
    """
    # Numbered 1, which a form may hold on its own and unconjugated in
    # complex mode.
    argument = ufl.Argument(space, 1, part)
    interpolation = ufl.Interpolate(argument, target)
    if processing.form:
        shape = interpolation.ufl_shape
        components = itertools.product(*(range(n) for n in shape))
        integrand = sum(interpolation[c] for c in components) if shape else interpolation
        # Restrictions are applied to an interpolation rather than into its
        # operand, and the scaling to the integrand, so a cell integral
        # processes it as an integral of any type does.
        form_data = _compute_form_data(
            integrand * ufl.Measure("dx", domain=target.ufl_domain()),
            processing.complex_mode,
            processing.do_cancel_jacobian_products,
        )
        expressions = [
            integral.integrand()
            for integral_data in form_data.integral_data
            for integral in integral_data.integrals
        ]
    else:
        expressions = [_process_expression(interpolation, processing.complex_mode)]
    (processed,) = {
        node
        for expression in expressions
        for node in unique_pre_traversal(expression)
        if isinstance(node, ufl.Interpolate)
    }
    return argument, renumber_indices(interpolated_expression(processed))


def _classify_interpolation(
    interpolation: ufl.Interpolate, processing: _Processing
) -> ufl.Argument | None:
    """Get the argument that a processed interpolation maps directly, if it maps one.

    After processing, the operand of an interpolation of an argument itself is
    the argument mapped by its own pull back and the inverse of the target
    element's, with the geometry lowered, so it is compared against an
    interpolation of the argument processed the same way. Cancelling Jacobian
    products in a form reduces it to the argument's reference value when the
    two pull backs are the same, which is compared against directly.

    Args:
        interpolation: A `ufl.Interpolate` from a processed expression.
        processing: How the expression has been processed.
    """
    argument = interpolated_argument(interpolation)
    if argument is not None:
        return argument
    arguments = interpolation_arguments(interpolation)
    if len(arguments) != 1:
        return None
    (argument,) = arguments
    source = argument.ufl_function_space()
    target = interpolation_target_space(interpolation)
    # The table of an argument itself maps between the reference elements,
    # which needs them on the same cells and with the same value shape. Only
    # then can the interpolation of the argument be built to compare against.
    if source.ufl_domain() != target.ufl_domain() or source.value_shape != target.value_shape:
        return None
    built, processed = _processed_argument(source, argument.part(), target, processing)
    operand = renumber_indices(interpolated_expression(interpolation))
    return argument if operand == ufl.replace(processed, {built: argument}) else None


def check_interpolation(interpolation: ufl.Interpolate, argument: ufl.Argument | None) -> None:
    """Check that an interpolation holding an argument can be turned into a table.

    Args:
        interpolation: A `ufl.Interpolate`.
        argument: The argument it maps directly, or None if it interpolates an
            expression merely linear in its arguments.
    """
    if argument is None:
        raise NotImplementedError(
            "Interpolating an expression that is linear in an argument is only "
            "supported when the expression is the argument itself."
        )
    # The argument's reference basis is pushed straight through the
    # interpolation, so the table is known at compile time. That needs the
    # Jacobian factors of the two pull backs to cancel and the block structures
    # to line up.
    target = interpolation_target_space(interpolation).ufl_element()
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


def interpolation_target_space(interpolation: ufl.Interpolate) -> ufl.FunctionSpace:
    """Get the function space that an interpolation lands in.

    Args:
        interpolation: A `ufl.Interpolate`.
    """
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
    expression: ufl.core.expr.Expr, scalar_type: npt.DTypeLike
) -> ufl.core.expr.Expr:
    """Analyzes and preprocesses expressions.

    Args:
        expression: The expression to process.
        scalar_type: The scalar type to compile for.
    """
    complex_mode = bool(np.issubdtype(scalar_type, np.complexfloating))
    expression = _process_expression(expression, complex_mode)
    replacer = IntermediateCoefficientReplacer(_Processing(form=False, complex_mode=complex_mode))
    return typing.cast(ufl.core.expr.Expr, replacer(expression))


def _process_expression(expression: ufl.core.expr.Expr, complex_mode: bool) -> ufl.core.expr.Expr:
    """Process an expression, as `ufl.algorithms.compute_form_data` does a form.

    Args:
        expression: The expression to process.
        complex_mode: Whether to process it for a complex scalar type.
    """
    preserve_geometry_types = (ufl.classes.Jacobian,)
    lowered = ufl.algorithms.apply_algebra_lowering.apply_algebra_lowering(expression)
    assert isinstance(lowered, ufl.core.expr.Expr)
    # An interpolation is linear in the expression it interpolates, so UFL
    # pushes a derivative into it.
    expression = apply_derivatives(lowered)
    # The same passes, in the same order, as `ufl.algorithms.compute_form_data`
    expression = ufl.algorithms.apply_function_pullbacks.apply_function_pullbacks(expression)
    expression = apply_geometry_lowering(expression, preserve_geometry_types)
    expression = apply_derivatives(expression)
    expression = apply_geometry_lowering(expression, preserve_geometry_types)
    expression = apply_derivatives(expression)

    # Remove complex nodes if scalar type is real valued
    if not complex_mode:
        expression = ufl.algorithms.remove_complex_nodes.remove_complex_nodes(expression)
    return expression


def _compute_form_data(
    form: ufl.Form, complex_mode: bool, do_cancel_jacobian_products: bool
) -> FormData:
    """Process a form with the options FFCx compiles forms with.

    Args:
        form: The form.
        complex_mode: Whether to process it for a complex scalar type.
        do_cancel_jacobian_products: Whether to cancel Jacobian products.
    """
    return ufl.algorithms.compute_form_data(
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
    complex_mode = bool(np.issubdtype(scalar_type, np.complexfloating))

    # Compute form metadata
    form_data = _compute_form_data(form, complex_mode, do_cancel_jacobian_products)

    # Shared between the integrals, which have been processed alike
    replacer = IntermediateCoefficientReplacer(
        _Processing(
            form=True,
            complex_mode=complex_mode,
            do_cancel_jacobian_products=do_cancel_jacobian_products,
        )
    )

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

            # Interpolations of an argument are classified. UFL has mapped
            # each onto the reference cell of its target element, and
            # represented what it reads in reference value.
            integrand = typing.cast(ufl.core.expr.Expr, replacer(integral.integrand()))

            # The elements of the arguments and coefficients the integrand
            # reads. UFL drops the coargument of an interpolation from an
            # expression, but keeps it for an integral, where its number
            # would clash with an argument's in a bilinear form.
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


class IntermediateCoefficientReplacer(DAGTraverser):
    """Classify the interpolations of an expression.

    An interpolation of an argument is left in place: it is a modified terminal,
    see `ffcx.ir.analysis.modified_terminals`, and its operand is rewritten into
    the argument's reference value, which `interpolated_argument` recognises
    from then on.
    """

    def __init__(self, processing: _Processing):
        """Initialise.

        Args:
            processing: How the expressions to replace in have been
                processed, to classify their interpolations by.
        """
        super().__init__()
        self._processing = processing

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
        if not interpolation_arguments(o):
            raise NotImplementedError(
                "Interpolating an expression that holds no argument is not supported."
            )
        argument = _classify_interpolation(o, self._processing)
        check_interpolation(o, argument)
        assert argument is not None
        # Interpolating the argument itself is a fixed map between the
        # reference elements, by the same pull back.
        return typing.cast(
            ufl.core.expr.Expr, o._ufl_expr_reconstruct_(ufl.classes.ReferenceValue(argument))
        )

    @process.register(ufl.core.expr.Expr)
    def _(self, o: ufl.core.expr.Expr) -> ufl.core.expr.Expr:
        """Handle anything else."""
        return self.reuse_if_untouched(o)
