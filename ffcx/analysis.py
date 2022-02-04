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
from collections import namedtuple

import numpy
import ufl

logger = logging.getLogger("ffcx")


ufl_data = namedtuple('ufl_data', ['form_data', 'unique_elements', 'element_numbers',
                                   'unique_coordinate_elements', 'expressions'])


def analyze_ufl_objects(ufl_objects: typing.List,
                        parameters: typing.Dict) -> ufl_data:
    """Analyze ufl object(s).

    Parameters
    ----------
    ufl_objects
    parameters
      FFCx parameters. These parameters take priority over all other set parameters.

    Returns
    -------
    form_datas
        Form_data objects
    unique_elements
        Unique elements across all forms
    element_numbers
        Mapping to unique numbers for all elements
    unique_coordinate_elements

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
            forms += [ufl_object]
        elif isinstance(ufl_object, ufl.FiniteElementBase):
            elements += [ufl_object]
        elif isinstance(ufl_object, ufl.Mesh):
            coordinate_elements += [ufl_object.ufl_coordinate_element()]
        elif isinstance(ufl_object[0], ufl.core.expr.Expr):
            original_expression = ufl_object[0]
            points = numpy.asarray(ufl_object[1])
            expressions += [(original_expression, points)]
        else:
            raise TypeError("UFL objects not recognised.")

    form_data = tuple(_analyze_form(form, parameters) for form in forms)
    for data in form_data:
        elements += data.unique_sub_elements
        coordinate_elements += data.coordinate_elements

    for original_expression, points in expressions:
        elements += list(ufl.algorithms.extract_elements(original_expression))
        processed_expression = _analyze_expression(original_expression, parameters)
        processed_expressions += [(processed_expression, points, original_expression)]

    elements += ufl.algorithms.analysis.extract_sub_elements(elements)

    # Sort elements so sub-elements come before mixed elements
    unique_elements = ufl.algorithms.sort_elements(set(elements))
    unique_coordinate_element_list = sorted(set(coordinate_elements), key=lambda x: repr(x))

    # Compute dict (map) from element to index
    element_numbers = {element: i for i, element in enumerate(unique_elements)}

    return ufl_data(form_data=form_data, unique_elements=unique_elements,
                    element_numbers=element_numbers,
                    unique_coordinate_elements=unique_coordinate_element_list,
                    expressions=processed_expressions)


def _analyze_expression(expression: ufl.core.expr.Expr, parameters: typing.Dict):
    """Analyzes and preprocesses expressions."""
    preserve_geometry_types = (ufl.classes.Jacobian, )
    expression = ufl.algorithms.apply_algebra_lowering.apply_algebra_lowering(expression)
    expression = ufl.algorithms.apply_derivatives.apply_derivatives(expression)
    expression = ufl.algorithms.apply_function_pullbacks.apply_function_pullbacks(expression)
    expression = ufl.algorithms.apply_geometry_lowering.apply_geometry_lowering(expression, preserve_geometry_types)
    expression = ufl.algorithms.apply_derivatives.apply_derivatives(expression)
    expression = ufl.algorithms.apply_geometry_lowering.apply_geometry_lowering(expression, preserve_geometry_types)
    expression = ufl.algorithms.apply_derivatives.apply_derivatives(expression)

    complex_mode = "_Complex" in parameters["scalar_type"]
    if not complex_mode:
        expression = ufl.algorithms.remove_complex_nodes.remove_complex_nodes(expression)

    return expression


def _analyze_form(form: ufl.form.Form, parameters: typing.Dict) -> ufl.algorithms.formdata.FormData:
    """Analyzes UFL form and attaches metadata.

    Parameters
    ----------
    form
    parameters

    Returns
    -------
    form_data -  Form data computed by UFL with metadata attached

    Note
    ----
    The main workload of this function is extraction of unique/default metadata
    from parameters, integral metadata or inherited from UFL
    (in case of quadrature degree)

    """
    if form.empty():
        raise RuntimeError(f"Form ({form}) seems to be zero: cannot compile it.")
    if _has_custom_integrals(form):
        raise RuntimeError(f"Form ({form}) contains unsupported custom integrals.")

    # Set default spacing for coordinate elements to be equispaced
    for n, i in enumerate(form._integrals):
        element = i._ufl_domain._ufl_coordinate_element
        if element._sub_element._variant is None and element.degree() > 2:
            sub_element = ufl.FiniteElement(
                element.family(), element.cell(), element.degree(), element.quadrature_scheme(),
                variant="equispaced")
            equi_element = ufl.VectorElement(sub_element)
            form._integrals[n]._ufl_domain._ufl_coordinate_element = equi_element

    # Check for complex mode
    complex_mode = "_Complex" in parameters["scalar_type"]

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

    # Determine unique quadrature degree, quadrature scheme and
    # precision per each integral data
    for id, integral_data in enumerate(form_data.integral_data):
        # Iterate through groups of integral data. There is one integral
        # data for all integrals with same domain, itype, subdomain_id
        # (but possibly different metadata).
        #
        # Quadrature degree and quadrature scheme must be the same for
        # all integrals in this integral data group, i.e. must be the
        # same for for the same (domain, itype, subdomain_id)

        # Extract precision
        p_default = -1
        precisions = set([integral.metadata().get("precision", p_default)
                          for integral in integral_data.integrals])
        precisions.discard(p_default)

        if len(precisions) == 1:
            p = precisions.pop()
        elif len(precisions) == 0:
            # Default precision
            p = numpy.finfo("double").precision + 1  # == 16
        else:
            raise RuntimeError("Only one precision allowed within integrals grouped by subdomain.")

        integral_data.metadata["precision"] = p

        qd_default = -1
        qr_default = "default"

        for i, integral in enumerate(integral_data.integrals):
            # Extract quadrature degree
            qd_metadata = integral.metadata().get("quadrature_degree", qd_default)
            pd_estimated = numpy.max(integral.metadata()["estimated_polynomial_degree"])
            if qd_metadata != qd_default:
                qd = qd_metadata
            else:
                qd = pd_estimated

            # Extract quadrature rule
            qr = integral.metadata().get("quadrature_rule", qr_default)

            logger.info(f"Integral {i}, integral group {id}:")
            logger.info(f"--- quadrature rule: {qr}")
            logger.info(f"--- quadrature degree: {qd}")
            logger.info(f"--- precision: {p}")

            # Update the old metadata
            metadata = integral.metadata()
            metadata.update({"quadrature_degree": qd, "quadrature_rule": qr, "precision": p})

            integral_data.integrals[i] = integral.reconstruct(metadata=metadata)

    return form_data


def _has_custom_integrals(o) -> bool:
    """Check for custom integrals."""
    if isinstance(o, ufl.integral.Integral):
        return o.integral_type() in ufl.custom_integral_types
    elif isinstance(o, ufl.classes.Form):
        return any(_has_custom_integrals(itg) for itg in o.integrals())
    elif isinstance(o, (list, tuple)):
        return any(_has_custom_integrals(itg) for itg in o)
    else:
        raise NotImplementedError
