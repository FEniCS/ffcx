# Copyright (C) 2007-2020 Anders Logg, Martin Alnaes, Kristian B. Oelgaard,
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

import logging
import typing

import basix.ufl
import numpy as np
import numpy.typing as npt
import ufl

logger = logging.getLogger("ffcx")


class UFLData(typing.NamedTuple):
    """UFL data."""

    # Tuple of ufl form data
    form_data: tuple[ufl.algorithms.formdata.FormData, ...]
    # List of unique elements
    unique_elements: list[basix.ufl._ElementBase]
    # Lookup table from each unique element to its index in `unique_elements`
    element_numbers: dict[basix.ufl._ElementBase, int]
    # List of unique coordinate elements
    unique_coordinate_elements: list[basix.ufl._ElementBase]
    # List of ufl Expressions as tuples (expression, points, original_expression)
    expressions: list[tuple[ufl.core.expr.Expr, npt.NDArray[np.float64], ufl.core.expr.Expr]]


def analyze_ufl_objects(
    ufl_objects: list[
        ufl.form.Form
        | ufl.AbstractFiniteElement
        | ufl.Mesh
        | tuple[ufl.core.expr.Expr, npt.NDArray[np.floating]]
    ],
    scalar_type: str,
) -> UFLData:
    """Analyze ufl object(s).

    Args:
        ufl_objects: UFL objects
        scalar_type: Scalar type that should be used for the analysis

    Returns:
        A data structure holding:
            form_datas: Form_data objects
            unique_elements: Unique elements across all forms and expressions
            element_numbers: Mapping to unique numbers for all elements
            unique_coordinate_elements: Unique coordinate elements across all forms and expressions
            expressions: List of all expressions after post-processing, with its evaluation points
                         and the original expression
    """
    logger.info(79 * "*")
    logger.info("Compiler stage 1: Analyzing UFL objects")
    logger.info(79 * "*")

    elements = []
    coordinate_elements = []

    # Group objects by types
    forms = []
    expressions = []
    processed_expressions = []

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


def _analyze_expression(expression: ufl.core.expr.Expr, scalar_type: str) -> ufl.core.expr.Expr:
    """Analyzes and preprocesses expressions."""
    preserve_geometry_types = (ufl.classes.Jacobian,)
    expression = ufl.algorithms.apply_algebra_lowering.apply_algebra_lowering(expression)
    expression = ufl.algorithms.apply_derivatives.apply_derivatives(expression)
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


def _analyze_form(form: ufl.form.Form, scalar_type: str) -> ufl.algorithms.formdata.FormData:
    """Analyzes UFL form and attaches metadata.

    Args:
        form: forms
        scalar_type: Scalar type used for form. This is used to simplify real valued forms

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
    for i in form._integrals:
        assert isinstance(i._ufl_domain._ufl_coordinate_element, basix.ufl._ElementBase)

    # Check for complex mode
    complex_mode = np.issubdtype(scalar_type, np.complexfloating)

    # Compute form metadata
    form_data = ufl.algorithms.compute_form_data(
        form,
        do_apply_function_pullbacks=True,
        do_apply_integral_scaling=True,
        do_apply_geometry_lowering=True,
        preserve_geometry_types=(ufl.classes.Jacobian,),
        do_apply_restrictions=True,
        do_append_everywhere_integrals=False,  # do not add dx integrals to dx(i) in UFL
        complex_mode=complex_mode,
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

        qd_default = -1
        qr_default = "default"

        for i, integral in enumerate(integral_data.integrals):
            metadata = integral.metadata()
            # If form contains a quadrature element, use the custom quadrature scheme
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
                qd_metadata = integral.metadata().get("quadrature_degree", qd_default)
                pd_estimated = np.max(integral.metadata()["estimated_polynomial_degree"])
                if qd_metadata != qd_default:
                    qd = qd_metadata
                else:
                    qd = pd_estimated

                # Extract quadrature rule
                qr = integral.metadata().get("quadrature_rule", qr_default)

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
    o: typing.Union[ufl.integral.Integral, ufl.classes.Form, list, tuple],
) -> bool:
    """Check for custom integrals."""
    if isinstance(o, ufl.integral.Integral):
        return o.integral_type() in ufl.custom_integral_types
    elif isinstance(o, ufl.classes.Form):
        return any(_has_custom_integrals(itg) for itg in o.integrals())
    elif isinstance(o, (list, tuple)):
        return any(_has_custom_integrals(itg) for itg in o)
    else:
        raise NotImplementedError
