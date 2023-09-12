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
from warnings import warn

import numpy as np
import numpy.typing as npt

import basix.ufl
import ufl
from ffcx.element_interface import convert_element

logger = logging.getLogger("ffcx")


class UFLData(typing.NamedTuple):
    form_data: typing.Tuple[ufl.algorithms.formdata.FormData, ...]  # Tuple of ufl form data
    unique_elements: typing.List[basix.ufl._ElementBase]  # List of unique elements
    # Lookup table from each unique element to its index in `unique_elements`
    element_numbers: typing.Dict[basix.ufl._ElementBase, int]
    unique_coordinate_elements: typing.List[basix.ufl._ElementBase]  # List of unique coordinate elements
    # List of ufl Expressions as tuples (expression, points, original_expression)
    expressions: typing.List[typing.Tuple[ufl.core.expr.Expr, npt.NDArray[np.float64], ufl.core.expr.Expr]]


def analyze_ufl_objects(ufl_objects: typing.List, options: typing.Dict) -> UFLData:
    """Analyze ufl object(s).

    Options
    ----------
    ufl_objects
    options
      FFCx options. These options take priority over all other set options.

    Returns a data structure holding
    -------
    form_datas
        Form_data objects
    unique_elements
        Unique elements across all forms and expressions
    element_numbers
        Mapping to unique numbers for all elements
    unique_coordinate_elements
        Unique coordinate elements across all forms and expressions
    expressions
        List of all expressions after post-processing, with its
        evaluation points and the original expression

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
        elif isinstance(ufl_object, ufl.FiniteElementBase):
            elements.append(convert_element(ufl_object))
        elif isinstance(ufl_object, ufl.Mesh):
            coordinate_elements.append(convert_element(ufl_object.ufl_coordinate_element()))
        elif isinstance(ufl_object[0], ufl.core.expr.Expr):
            original_expression = ufl_object[0]
            points = np.asarray(ufl_object[1])
            expressions.append((original_expression, points))
        else:
            raise TypeError("UFL objects not recognised.")

    form_data = tuple(_analyze_form(form, options) for form in forms)
    for data in form_data:
        elements += [convert_element(e) for e in data.unique_sub_elements]
        coordinate_elements += [convert_element(e) for e in data.coordinate_elements]

    for original_expression, points in expressions:
        elements += [convert_element(e) for e in ufl.algorithms.extract_elements(original_expression)]
        processed_expression = _analyze_expression(original_expression, options)
        processed_expressions += [(processed_expression, points, original_expression)]

    elements += ufl.algorithms.analysis.extract_sub_elements(elements)

    # Sort elements so sub-elements come before mixed elements
    unique_elements = ufl.algorithms.sort_elements(set(elements))
    unique_coordinate_element_list = sorted(set(coordinate_elements), key=lambda x: repr(x))

    for e in unique_elements:
        assert isinstance(e, basix.ufl._ElementBase)

    # Compute dict (map) from element to index
    element_numbers = {element: i for i, element in enumerate(unique_elements)}

    return UFLData(form_data=form_data, unique_elements=unique_elements, element_numbers=element_numbers,
                   unique_coordinate_elements=unique_coordinate_element_list, expressions=processed_expressions)


def _analyze_expression(expression: ufl.core.expr.Expr, options: typing.Dict):
    """Analyzes and preprocesses expressions."""
    preserve_geometry_types = (ufl.classes.Jacobian, )
    expression = ufl.algorithms.apply_algebra_lowering.apply_algebra_lowering(expression)
    expression = ufl.algorithms.apply_derivatives.apply_derivatives(expression)
    expression = ufl.algorithms.apply_function_pullbacks.apply_function_pullbacks(expression)
    expression = ufl.algorithms.apply_geometry_lowering.apply_geometry_lowering(expression, preserve_geometry_types)
    expression = ufl.algorithms.apply_derivatives.apply_derivatives(expression)
    expression = ufl.algorithms.apply_geometry_lowering.apply_geometry_lowering(expression, preserve_geometry_types)
    expression = ufl.algorithms.apply_derivatives.apply_derivatives(expression)

    complex_mode = "_Complex" in options["scalar_type"]
    if not complex_mode:
        expression = ufl.algorithms.remove_complex_nodes.remove_complex_nodes(expression)

    return expression


def _analyze_form(form: ufl.form.Form, options: typing.Dict) -> ufl.algorithms.formdata.FormData:
    """Analyzes UFL form and attaches metadata.

    Args:
        form: forms
        options: options

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

    # Set default spacing for coordinate elements to be equispaced
    for n, i in enumerate(form._integrals):
        element = i._ufl_domain._ufl_coordinate_element
        if not isinstance(element, basix.ufl._ElementBase) and element.degree() > 2:
            warn("UFL coordinate elements using elements not created via Basix may not work with DOLFINx")

    # Check for complex mode
    complex_mode = "_Complex" in options["scalar_type"]

    # Compute form metadata
    form_data = ufl.algorithms.compute_form_data(
        form,
        do_apply_function_pullbacks=True,
        do_apply_integral_scaling=True,
        do_apply_geometry_lowering=True,
        preserve_geometry_types=(ufl.classes.Jacobian,),
        do_apply_restrictions=True,
        do_append_everywhere_integrals=False,  # do not add dx integrals to dx(i) in UFL
        complex_mode=complex_mode)

    # If form contains a quadrature element, use the custom quadrature scheme
    custom_q = None
    for e in form_data.unique_elements:
        e = convert_element(e)
        if e.has_custom_quadrature:
            if custom_q is None:
                custom_q = e.custom_quadrature()
            else:
                assert np.allclose(e._points, custom_q[0])
                assert np.allclose(e._weights, custom_q[1])

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
                metadata.update({"quadrature_points": custom_q[0], "quadrature_weights": custom_q[1],
                                 "quadrature_rule": "custom"})

            integral_data.integrals[i] = integral.reconstruct(metadata=metadata)

    return form_data


def _has_custom_integrals(o: typing.Union[ufl.integral.Integral, ufl.classes.Form, list, tuple]) -> bool:
    """Check for custom integrals."""
    if isinstance(o, ufl.integral.Integral):
        return o.integral_type() in ufl.custom_integral_types
    elif isinstance(o, ufl.classes.Form):
        return any(_has_custom_integrals(itg) for itg in o.integrals())
    elif isinstance(o, (list, tuple)):
        return any(_has_custom_integrals(itg) for itg in o)
    else:
        raise NotImplementedError
