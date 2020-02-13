# Copyright (C) 2007-2020 Anders Logg, Martin Alnaes, Kristian B. Oelgaard,
#                         Michal Habera and others
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Compiler stage 1: Analysis.

This module implements the analysis/preprocessing of variational forms,
including automatic selection of elements, degrees and form
representation type.
"""

import logging
import os
import typing
import warnings
from collections import namedtuple

import numpy

import ufl
from ffcx.parameters import FFCX_PARAMETERS
logger = logging.getLogger(__name__)


ufl_data = namedtuple('ufl_data', ['form_data', 'unique_elements', 'element_numbers',
                                   'unique_coordinate_elements', 'expressions'])


def analyze_ufl_objects(ufl_objects: typing.Union[typing.List[ufl.form.Form], typing.List[ufl.FiniteElement],
                                                  typing.List],
                        parameters: typing.Dict) -> ufl_data:
    """Analyze ufl object(s).

    Parameters
    ----------
    ufl_objects
    parameters

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
    logger.info("Compiler stage 1: Analyzing UFL objects")

    form_data = ()
    unique_elements = set()
    unique_coordinate_elements = set()
    expressions = []

    # FIXME: This assumes that forms come before elements in
    # ufl_objects? Is this reasonable?
    if isinstance(ufl_objects[0], ufl.form.Form):
        forms = ufl_objects
        form_data = tuple(_analyze_form(form, parameters) for form in forms)

        # Extract unique elements across forms
        for data in form_data:
            unique_elements.update(data.unique_sub_elements)

        # Extract uniquecoordinate elements across forms
        for data in form_data:
            unique_coordinate_elements.update(data.coordinate_elements)
    elif isinstance(ufl_objects[0], ufl.FiniteElementBase):
        # Extract unique (sub)elements
        elements = ufl_objects
        unique_elements.update(ufl.algorithms.analysis.extract_sub_elements(elements))
    elif isinstance(ufl_objects[0], ufl.Mesh):
        # Extract unique (sub)elements
        meshes = ufl_objects
        unique_coordinate_elements = [mesh.ufl_coordinate_element() for mesh in meshes]
    elif isinstance(ufl_objects[0], tuple) and isinstance(ufl_objects[0][0], ufl.core.expr.Expr):
        for expression in ufl_objects:
            original_expression = expression[0]
            points = expression[1]
            expression = expression[0]

            unique_elements.update(ufl.algorithms.extract_elements(expression))
            unique_elements.update(ufl.algorithms.extract_sub_elements(unique_elements))

            expression = _analyze_expression(expression, parameters)
            expressions.append((expression, points, original_expression))
    else:
        raise TypeError("UFL objects not recognised.")

    # Make sure coordinate elements and their subelements are included
    unique_elements.update(ufl.algorithms.analysis.extract_sub_elements(unique_coordinate_elements))

    # Sort elements so sub-elements come before mixed elements
    unique_elements = ufl.algorithms.sort_elements(unique_elements)
    unique_coordinate_elements = sorted(unique_coordinate_elements, key=lambda x: repr(x))

    # Compute dict (map) from element to index
    element_numbers = {element: i for i, element in enumerate(unique_elements)}

    return ufl_data(form_data=form_data, unique_elements=unique_elements,
                    element_numbers=element_numbers,
                    unique_coordinate_elements=unique_coordinate_elements,
                    expressions=expressions)


def _analyze_expression(expression: ufl.core.expr.Expr, parameters: typing.Dict):
    """Analyzes and preprocesses expressions."""
    preserve_geometry_types = (ufl.CellVolume, ufl.FacetArea)
    expression = ufl.algorithms.apply_algebra_lowering.apply_algebra_lowering(expression)
    expression = ufl.algorithms.apply_derivatives.apply_derivatives(expression)
    expression = ufl.algorithms.apply_function_pullbacks.apply_function_pullbacks(expression)
    expression = ufl.algorithms.apply_geometry_lowering.apply_geometry_lowering(expression, preserve_geometry_types)
    expression = ufl.algorithms.apply_derivatives.apply_derivatives(expression)
    expression = ufl.algorithms.apply_geometry_lowering.apply_geometry_lowering(expression, preserve_geometry_types)
    expression = ufl.algorithms.apply_derivatives.apply_derivatives(expression)

    complex_mode = "complex" in parameters.get("scalar_type", "double")
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
        raise RuntimeError("Form ({}) seems to be zero: cannot compile it.".format(str(form)))
    if _has_custom_integrals(form):
        raise RuntimeError("Form ({}) contains unsupported custom integrals.".format(str(form)))

    # ---- Extract representation across all integrals in this form
    #
    # The priority of representation determination is following
    #
    # 1. Environment variable FFCX_FORCE_REPRESENTATION
    # 2. parameters["representation"]
    # 3. specified in metadata of integral
    representations = set(
        integral.metadata().get("representation", "auto") for integral in form.integrals())

    # Remove "auto" to see representations set by user
    if parameters["representation"] in ("uflacs", "tsfc"):
        representation = parameters["representation"]
    elif len(representations - {"auto"}) == 1:
        # User has set just one
        representation = (representations - {"auto"}).pop()
    elif representations == {"auto"}:
        # If user didn't set any default to uflacs
        representation = "uflacs"
    else:
        raise RuntimeError("Cannot mix uflacs and tsfc representations in a single form.")

    # Override representation with environment variable
    forced_r = os.environ.get("FFCX_FORCE_REPRESENTATION")
    if forced_r:
        warnings.warn("Representation: forced by $FFCX_FORCE_REPRESENTATION to '{}'".format(forced_r))
        representation = forced_r

    logger.info("Determined representation '{}' for form {}.".format(representation, str(form)))

    # Check for complex mode
    complex_mode = "complex" in parameters.get("scalar_type", "double")

    # Compute form metadata
    if representation == "uflacs":
        form_data = ufl.algorithms.compute_form_data(
            form,
            do_apply_function_pullbacks=True,
            do_apply_integral_scaling=True,
            do_apply_geometry_lowering=True,
            preserve_geometry_types=(ufl.classes.Jacobian, ),
            do_apply_restrictions=True,
            do_append_everywhere_integrals=False,  # do not add dx integrals to dx(i) in UFL
            complex_mode=complex_mode)
    elif representation == "tsfc":
        # TSFC provides compute_form_data wrapper using correct kwargs
        from tsfc.ufl_utils import compute_form_data as tsfc_compute_form_data
        form_data = tsfc_compute_form_data(form, complex_mode=complex_mode)
    else:
        raise RuntimeError("Unexpected representation \"{}\" for form preprocessing.".format(representation))

    # Attach common representation to FormData. Common to all integrals
    # in form
    form_data.representation = representation

    # Determine unique quadrature degree, quadrature scheme and
    # precision per each integral data
    for integral_data in form_data.integral_data:
        # Iterate through groups of integral data. There is one integral
        # data for all integrals with same domain, itype, subdomain_id
        # (but possibly different metadata).
        #
        # Quadrature degree and quadrature scheme must be the same for
        # all integrals in this integral data group, i.e. must be the
        # same for for the same (domain, itype, subdomain_id)

        # ----- Extract common quadrature degree
        #
        # The priority of quadrature degree determination is following
        #
        # 1. Parameters["quadrature_degree"]
        # 2. Specified in metadata of integral
        # 3. Estimated by UFL
        default = FFCX_PARAMETERS["quadrature_degree"]
        quadrature_degrees = set([integral.metadata().get("quadrature_degree", default)
                                  for integral in integral_data.integrals])
        quadrature_degrees.discard(default)

        if parameters["quadrature_degree"] != default:
            # Quadrature degree is forced by FFCX parameters
            qd = parameters["quadrature_degree"]
        elif len(quadrature_degrees) == 1:
            qd = quadrature_degrees.pop()
        elif len(quadrature_degrees) == 0:
            # If there are more integrals in this integral data group
            # and UFL estimated different degrees we pick maximum

            # Find all estimated polynomial degrees by UFL for all integrals in this
            # integral data group
            estimated_quadrature_degrees = [integral.metadata()["estimated_polynomial_degree"]
                                            for integral in integral_data.integrals]
            qd = max(estimated_quadrature_degrees)

            # The quadrature degree from UFL can be very high for some
            # integrals.  Print warning if number of quadrature points
            # exceeds 100.
            tdim = integral_data.domain.topological_dimension()
            num_points = ((qd + 1 + 1) // 2)**tdim
            if num_points >= 100:
                warnings.warn(
                    "Number of integration points per cell is: {}. Consider using 'quadrature_degree' "
                    "to reduce number.".format(num_points))
        elif len(quadrature_degrees) > 1:
            raise RuntimeError("Only one quadrature degree allowed within integrals grouped by subdomain.")
        else:
            raise RuntimeError("Unable to determine quadrature degree.")

        # ----- Extract common quadrature rule
        #
        # The priority of quadrature rule determination is following
        #
        # 1. parameters["quadrature_rule"]
        # 2. specified in metadata of integral
        default = FFCX_PARAMETERS["quadrature_rule"]
        quadrature_rules = set([integral.metadata().get("quadrature_rule", default)
                                for integral in integral_data.integrals])
        quadrature_rules.discard(default)

        if parameters["quadrature_rule"] != default:
            qr = parameters["quadrature_rule"]
        elif len(quadrature_rules) == 1:
            qr = quadrature_rules.pop()
        elif len(quadrature_rules) == 0:
            qr = "default"
        elif len(quadrature_rules) > 1:
            raise RuntimeError("Only one quadrature scheme allowed within integrals grouped by subdomain.")
        else:
            raise RuntimeError("Unable to determine quadrature rule.")

        # ----- Extract precision
        #
        # The priority of precision determination is following
        #
        # 1. parameters["precision"]
        # 2. specified in metadata of integral
        default = FFCX_PARAMETERS["precision"]
        precisions = set([integral.metadata().get("precision", default)
                          for integral in integral_data.integrals])
        precisions.discard(default)

        if parameters["precision"] != default:
            p = parameters["precision"]
        elif len(precisions) == 1:
            p = precisions.pop()
        elif len(precisions) == 0:
            p = numpy.finfo("double").precision + 1  # == 16
        elif len(precisions) > 1:
            raise RuntimeError("Only one precision allowed within integrals grouped by subdomain.")
        else:
            raise RuntimeError("Unable to determine quadrature degree.")

        integral_data.metadata["quadrature_degree"] = qd
        integral_data.metadata["quadrature_rule"] = qr
        integral_data.metadata["precision"] = p

        # Reconstruct integrals to avoid modifying the input integral,
        # which would affect the signature computation if the integral
        # was used again in the user program.  Modifying attributes of
        # form_data.integral_data is less problematic since it's
        # lifetime is internal to the form compiler pipeline.
        for i, integral in enumerate(integral_data.integrals):
            integral_data.integrals[i] = integral.reconstruct(
                metadata={"quadrature_degree": qd, "quadrature_rule": qr, "precision": p})

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
