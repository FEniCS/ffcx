# Copyright (C) 2013-2017 Martin Sandve Alnæs
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import logging
import hashlib

import numpy

import ufl
from ffcx.fiatinterface import create_element
from ffcx.ir.representationutils import create_quadrature_points_and_weights
from ffcx.ir.uflacs.build_uflacs_ir import build_uflacs_ir
from ufl.sorting import sorted_expr_sum
from ufl.classes import Integral


logger = logging.getLogger(__name__)


class QuadratureRule:
    def __init__(self, points, weights):
        self.points = points
        self.weights = weights
        self._hash = None

    def __hash__(self):
        if self._hash is None:
            self.hash_obj = hashlib.sha1(self.points)
            self._hash = int(self.hash_obj.hexdigest(), 32)
        return self._hash

    def __eq__(self, other):
        return numpy.allclose(self.points, other.points) and numpy.allclose(self.weights, other.weights)

    def id(self):
        return self.hash_obj.hexdigest()[-3:]


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

    weights = numpy.array([1.0] * points.shape[0])
    rule = QuadratureRule(points, weights)
    integrands = {rule: expression}

    uflacs_ir = build_uflacs_ir(cell, ir["integral_type"], ir["entitytype"], integrands, tensor_shape,
                                parameters, visualise)

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

    # Group integrands with the same quadrature rule
    grouped_integrands = {}
    for integral in itg_data.integrals:
        md = integral.metadata() or {}
        scheme = md["quadrature_rule"]
        degree = md["quadrature_degree"]

        (points, weights) = create_quadrature_points_and_weights(integral_type, cell, degree,
                                                                 scheme)
        points = numpy.asarray(points)
        weights = numpy.asarray(weights)

        rule = QuadratureRule(points, weights)

        if rule not in grouped_integrands:
            grouped_integrands[rule] = []

        grouped_integrands[rule].append(integral.integrand())

    sorted_integrals = {}
    for rule, integrands in grouped_integrands.items():
        integrands_summed = sorted_expr_sum(integrands)

        integral_new = Integral(integrands_summed, itg_data.integral_type, itg_data.domain,
                                itg_data.subdomain_id, {}, None)
        sorted_integrals[rule] = integral_new

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
    integrands = {rule: integral.integrand() for rule, integral in sorted_integrals.items()}

    # Build the more uflacs-specific intermediate representation
    uflacs_ir = build_uflacs_ir(itg_data.domain.ufl_cell(), itg_data.integral_type,
                                ir["entitytype"], integrands, ir["tensor_shape"],
                                parameters, visualise)

    ir.update(uflacs_ir)

    return ir
