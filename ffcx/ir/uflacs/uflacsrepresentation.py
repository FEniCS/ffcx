# Copyright (C) 2013-2017 Martin Sandve Alnæs
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import logging

import numpy

import ufl
from ffcx.fiatinterface import create_element
from ffcx.ir.uflacs.build_uflacs_ir import build_uflacs_ir
from ffcx.ir.uflacs.tools import accumulate_integrals, compute_quadrature_rules

logger = logging.getLogger(__name__)


def compute_expression_ir(expression, analysis, parameters, visualise):
    """Compute IR for expression.

    Parameters
    ----------
    expression
        Triple of (UFL expression, array of evaluation points, original UFL expression).

    Note
    ----
    Original UFL expression is needed to compute original positions of coefficients.

    """
    logger.info("Computing uflacs representation of expression")

    original_expression = expression[2]
    points = expression[1]
    expression = expression[0]

    num_points = points.shape[0]
    weights = numpy.array([1.0] * num_points)

    cell = expression.ufl_domain().ufl_cell()

    ir = {}

    # Prepare dimensions of all unique element in expression, including
    # elements for arguments, coefficients and coordinate mappings
    ir["element_dimensions"] = {
        ufl_element: create_element(ufl_element).space_dimension()
        for ufl_element in analysis.unique_elements
    }

    # Extract dimensions for elements of arguments only
    arguments = ufl.algorithms.extract_arguments(expression)
    argument_elements = tuple(f.ufl_element() for f in arguments)
    argument_dimensions = [
        ir["element_dimensions"][ufl_element] for ufl_element in argument_elements
    ]

    tensor_shape = argument_dimensions
    ir["tensor_shape"] = tensor_shape

    ir["expression_shape"] = list(expression.ufl_shape)

    coefficients = ufl.algorithms.extract_coefficients(expression)
    coefficient_numbering = {}
    for i, coeff in enumerate(coefficients):
        coefficient_numbering[coeff] = i

    # Add coefficient numbering to IR
    ir["coefficient_numbering"] = coefficient_numbering

    original_coefficient_positions = []
    original_coefficients = ufl.algorithms.extract_coefficients(original_expression)
    for coeff in coefficients:
        original_coefficient_positions.append(original_coefficients.index(coeff))

    ir["original_coefficient_positions"] = original_coefficient_positions

    coefficient_elements = tuple(f.ufl_element() for f in coefficients)

    offsets = {}
    _offset = 0
    for i, el in enumerate(coefficient_elements):
        offsets[coefficients[i]] = _offset
        _offset += ir["element_dimensions"][el]

    # Copy offsets also into IR
    ir["coefficient_offsets"] = offsets

    ir["integral_type"] = "expression"
    ir["entitytype"] = "cell"

    # Build offsets for Constants
    original_constant_offsets = {}
    _offset = 0
    for constant in ufl.algorithms.analysis.extract_constants(expression):
        original_constant_offsets[constant] = _offset
        _offset += numpy.product(constant.ufl_shape, dtype=numpy.int)

    ir["original_constant_offsets"] = original_constant_offsets

    ir["points"] = points

    integrands = {num_points: expression}
    quadrature_rules = {num_points: (points, weights)}

    uflacs_ir = build_uflacs_ir(cell, ir["integral_type"], ir["entitytype"], integrands, tensor_shape,
                                quadrature_rules, parameters, visualise)

    ir.update(uflacs_ir)

    return ir


def compute_integral_ir(ir: dict,
                        itg_data: ufl.algorithms.domain_analysis.IntegralData,
                        form_data: ufl.algorithms.formdata.FormData,
                        element_numbers: dict, parameters: dict,
                        visualise: bool):
    """Compute intermediate represention of integral."""

    logger.info("Computing uflacs representation for integral")

    # Get element space dimensions
    unique_elements = element_numbers.keys()
    ir["element_dimensions"] = {
        ufl_element: create_element(ufl_element).space_dimension()
        for ufl_element in unique_elements
    }

    ir["element_ids"] = {
        ufl_element: i
        for i, ufl_element in enumerate(unique_elements)
    }

    # Create dimensions of primary indices, needed to reset the argument
    # 'A' given to tabulate_tensor() by the assembler.
    argument_dimensions = [
        ir["element_dimensions"][ufl_element] for ufl_element in form_data.argument_elements
    ]

    # Compute shape of element tensor
    if ir["integral_type"] == "interior_facet":
        ir["tensor_shape"] = [2 * dim for dim in argument_dimensions]
    else:
        ir["tensor_shape"] = argument_dimensions

    integral_type = itg_data.integral_type
    cell = itg_data.domain.ufl_cell()

    # Collect the quadrature rules occur in integrals
    rules = set()
    for integral in itg_data.integrals:
        md = integral.metadata() or {}
        scheme = md["quadrature_rule"]
        degree = md["quadrature_degree"]
        rules.add((scheme, degree))
    quadrature_integral_type = integral_type

    # Compute actual points and weights
    quadrature_rules, quadrature_rule_sizes = compute_quadrature_rules(
        rules, quadrature_integral_type, cell)

    # Store quadrature rules in format {num_points: (points, weights)}.
    # There is an assumption that schemes with the same number of points
    # are identical.
    ir["quadrature_rules"] = quadrature_rules

    # Group and accumulate integrals on the format {num_points: integral
    # data}
    sorted_integrals = accumulate_integrals(itg_data, quadrature_rule_sizes)

    # TODO: See if coefficient_numbering can be removed
    # Build coefficient numbering for UFC interface here, to avoid
    # renumbering in UFL and application of replace mapping
    coefficient_numbering = {}
    for i, f in enumerate(form_data.reduced_coefficients):
        coefficient_numbering[f] = i

    # Add coefficient numbering to IR
    ir["coefficient_numbering"] = coefficient_numbering

    index_to_coeff = sorted([(v, k) for k, v in coefficient_numbering.items()])
    offsets = {}
    width = 2 if integral_type in ("interior_facet") else 1
    _offset = 0
    for k, el in zip(index_to_coeff, form_data.coefficient_elements):
        offsets[k[1]] = _offset
        _offset += width * ir["element_dimensions"][el]

    # Copy offsets also into IR
    ir["coefficient_offsets"] = offsets

    # Build offsets for Constants
    original_constant_offsets = {}
    _offset = 0
    for constant in form_data.original_form.constants():
        original_constant_offsets[constant] = _offset
        _offset += numpy.product(constant.ufl_shape, dtype=numpy.int)

    ir["original_constant_offsets"] = original_constant_offsets

    ir["precision"] = itg_data.metadata["precision"]

    # Create map from number of quadrature points -> integrand
    integrands = {num_points: integral.integrand() for num_points, integral in sorted_integrals.items()}

    # Build the more uflacs-specific intermediate representation
    uflacs_ir = build_uflacs_ir(itg_data.domain.ufl_cell(), itg_data.integral_type,
                                ir["entitytype"], integrands, ir["tensor_shape"],
                                quadrature_rules, parameters, visualise)

    ir.update(uflacs_ir)

    return ir
