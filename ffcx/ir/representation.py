# Copyright (C) 2009-2020 Anders Logg, Martin Sandve Alnæs, Marie E. Rognes,
# Kristian B. Oelgaard, Matthew W. Scroggs, Chris Richardson, and others
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Compiler stage 2: Code representation.

Module computes intermediate representations of forms, elements and
dofmaps. For each UFC function, we extract the data needed for code
generation at a later stage.

The representation should conform strictly to the naming and order of
functions in UFC. Thus, for code generation of the function "foo", one
should only need to use the data stored in the intermediate
representation under the key "foo".
"""

import itertools
import logging
import warnings
from collections import namedtuple

import numpy
import ufl
from ffcx import naming
from ffcx.element_interface import create_element
from ffcx.ir.integral import compute_integral_ir
from ffcx.ir.representationutils import (QuadratureRule,
                                         create_quadrature_points_and_weights)
from ufl.classes import Integral
from ufl.sorting import sorted_expr_sum

logger = logging.getLogger("ffcx")

ir_form = namedtuple('ir_form', [
    'id', 'name', 'signature', 'rank', 'num_coefficients', 'num_constants',
    'name_from_uflfile', 'function_spaces', 'original_coefficient_position',
    'coefficient_names', 'constant_names', 'finite_elements',
    'dofmaps', 'integral_names', 'subdomain_ids'])
ir_element = namedtuple('ir_element', [
    'id', 'name', 'signature', 'cell_shape', 'topological_dimension',
    'geometric_dimension', 'space_dimension', 'value_shape', 'reference_value_shape', 'degree',
    'family', 'num_sub_elements', 'block_size', 'sub_elements', 'element_type', 'entity_dofs',
    'lagrange_variant', 'dpc_variant', 'basix_family', 'basix_cell', 'discontinuous', 'custom_element'])
ir_dofmap = namedtuple('ir_dofmap', [
    'id', 'name', 'signature', 'num_global_support_dofs', 'num_element_support_dofs', 'num_entity_dofs',
    'tabulate_entity_dofs', 'num_entity_closure_dofs', 'tabulate_entity_closure_dofs', 'num_sub_dofmaps',
    'sub_dofmaps', 'block_size'])
ir_integral = namedtuple('ir_integral', [
    'integral_type', 'subdomain_id', 'rank', 'geometric_dimension', 'topological_dimension', 'entitytype',
    'num_facets', 'num_vertices', 'enabled_coefficients', 'element_dimensions',
    'element_ids', 'tensor_shape', 'coefficient_numbering', 'coefficient_offsets',
    'original_constant_offsets', 'params', 'cell_shape', 'unique_tables', 'unique_table_types',
    'table_dofmaps', 'integrand', 'name', 'precision', 'needs_facet_permutations', 'coordinate_element'])
ir_expression = namedtuple('ir_expression', [
    'name', 'element_dimensions', 'params', 'unique_tables', 'unique_table_types', 'integrand',
    'table_dofmaps', 'coefficient_numbering', 'coefficient_offsets',
    'integral_type', 'entitytype', 'tensor_shape', 'expression_shape', 'original_constant_offsets',
    'original_coefficient_positions', 'points', 'coefficient_names', 'constant_names', 'needs_facet_permutations',
    'function_spaces', 'name_from_uflfile'])
ir_custom_element = namedtuple('ir_custom_element', [
    'cell_type', 'value_shape', 'wcoeffs', 'x', 'M', 'map_type',
    'discontinuous', 'highest_complete_degree', 'highest_degree'])

ir_data = namedtuple('ir_data', ['elements', 'dofmaps', 'integrals', 'forms', 'expressions'])


def compute_ir(analysis, object_names, prefix, parameters, visualise):
    """Compute intermediate representation."""
    logger.info(79 * "*")
    logger.info("Compiler stage 2: Computing intermediate representation of objects")
    logger.info(79 * "*")

    # Compute object names
    # NOTE: This is done here for performance reasons, because repeated calls
    # within each IR computation would be expensive due to UFL signature computations
    finite_element_names = {e: naming.finite_element_name(e, prefix) for e in analysis.unique_elements}
    dofmap_names = {e: naming.dofmap_name(e, prefix) for e in analysis.unique_elements}
    integral_names = {}
    form_names = {}
    for fd_index, fd in enumerate(analysis.form_data):
        form_names[fd_index] = naming.form_name(fd.original_form, fd_index, prefix)
        for itg_index, itg_data in enumerate(fd.integral_data):
            integral_names[(fd_index, itg_index)] = naming.integral_name(fd.original_form, itg_data.integral_type,
                                                                         fd_index, itg_data.subdomain_id, prefix)

    ir_elements = [
        _compute_element_ir(e, analysis.element_numbers, finite_element_names)
        for e in analysis.unique_elements
    ]

    ir_dofmaps = [
        _compute_dofmap_ir(e, analysis.element_numbers, dofmap_names)
        for e in analysis.unique_elements
    ]

    irs = [
        _compute_integral_ir(fd, i, analysis.element_numbers, integral_names, finite_element_names,
                             parameters, visualise)
        for (i, fd) in enumerate(analysis.form_data)
    ]
    ir_integrals = list(itertools.chain(*irs))

    ir_forms = [
        _compute_form_ir(fd, i, prefix, form_names, integral_names, analysis.element_numbers, finite_element_names,
                         dofmap_names, object_names)
        for (i, fd) in enumerate(analysis.form_data)
    ]

    ir_expressions = [_compute_expression_ir(expr, i, prefix, analysis, parameters, visualise, object_names,
                                             finite_element_names, dofmap_names)
                      for i, expr in enumerate(analysis.expressions)]

    return ir_data(elements=ir_elements, dofmaps=ir_dofmaps,
                   integrals=ir_integrals, forms=ir_forms,
                   expressions=ir_expressions)


def _compute_element_ir(ufl_element, element_numbers, finite_element_names):
    """Compute intermediate representation of element."""
    logger.info(f"Computing IR for element {ufl_element}")

    # Create basix elements
    basix_element = create_element(ufl_element)
    cell = ufl_element.cell()
    cellname = cell.cellname()

    # Store id
    ir = {"id": element_numbers[ufl_element]}
    ir["name"] = finite_element_names[ufl_element]

    # Compute data for each function
    ir["signature"] = repr(ufl_element)
    ir["cell_shape"] = cellname
    ir["topological_dimension"] = cell.topological_dimension()
    ir["geometric_dimension"] = cell.geometric_dimension()
    ir["space_dimension"] = basix_element.dim
    ir["element_type"] = basix_element.element_type
    ir["lagrange_variant"] = basix_element.lagrange_variant
    ir["dpc_variant"] = basix_element.dpc_variant
    ir["basix_family"] = basix_element.element_family
    ir["basix_cell"] = basix_element.cell_type
    ir["discontinuous"] = basix_element.discontinuous
    ir["degree"] = ufl_element.degree()
    ir["family"] = ufl_element.family()
    ir["value_shape"] = ufl_element.value_shape()
    ir["reference_value_shape"] = ufl_element.reference_value_shape()

    ir["num_sub_elements"] = ufl_element.num_sub_elements()
    ir["sub_elements"] = [finite_element_names[e] for e in ufl_element.sub_elements()]

    if hasattr(basix_element, "block_size"):
        ir["block_size"] = basix_element.block_size
        ufl_element = ufl_element.sub_elements()[0]
        basix_element = create_element(ufl_element)
    else:
        ir["block_size"] = 1

    ir["entity_dofs"] = basix_element.entity_dofs

    if basix_element.is_custom_element:
        ir["custom_element"] = _compute_custom_element_ir(basix_element.element)
    else:
        ir["custom_element"] = None

    return ir_element(**ir)


def _compute_custom_element_ir(basix_element):
    """Compute intermediate representation of a custom Basix element."""
    ir = {}
    ir["cell_type"] = basix_element.cell_type
    ir["value_shape"] = basix_element.value_shape
    ir["wcoeffs"] = basix_element.wcoeffs
    ir["x"] = basix_element.x
    ir["M"] = basix_element.M
    ir["map_type"] = basix_element.map_type
    ir["discontinuous"] = basix_element.discontinuous
    ir["highest_complete_degree"] = basix_element.highest_complete_degree
    ir["highest_degree"] = basix_element.highest_degree

    return ir_custom_element(**ir)


def _compute_dofmap_ir(ufl_element, element_numbers, dofmap_names):
    """Compute intermediate representation of dofmap."""
    logger.info(f"Computing IR for dofmap of {ufl_element}")

    # Create basix elements
    basix_element = create_element(ufl_element)

    # Store id
    ir = {"id": element_numbers[ufl_element]}
    ir["name"] = dofmap_names[ufl_element]

    # Compute data for each function
    ir["signature"] = "FFCx dofmap for " + repr(ufl_element)
    ir["sub_dofmaps"] = [dofmap_names[e] for e in ufl_element.sub_elements()]
    ir["num_sub_dofmaps"] = ufl_element.num_sub_elements()

    if hasattr(basix_element, "block_size"):
        ir["block_size"] = basix_element.block_size
        basix_element = basix_element.sub_element
    else:
        ir["block_size"] = 1

    # Precompute repeatedly used items
    for i in basix_element.num_entity_dofs:
        # FIXME: this assumes the same number of DOFs on each entity of the same dim: this
        # assumption will not be true for prisms and pyramids
        if max(i) != min(i):
            raise RuntimeError("Elements with different numbers of DOFs on subentities of the same dimension"
                               " are not yet supported in FFCx.")

    num_dofs_per_entity = [i[0] for i in basix_element.num_entity_dofs]
    ir["num_entity_dofs"] = num_dofs_per_entity
    ir["tabulate_entity_dofs"] = (basix_element.entity_dofs, num_dofs_per_entity)

    num_dofs_per_entity_closure = [i[0] for i in basix_element.num_entity_closure_dofs]
    ir["num_entity_closure_dofs"] = num_dofs_per_entity_closure
    ir["tabulate_entity_closure_dofs"] = (basix_element.entity_closure_dofs, num_dofs_per_entity_closure)

    ir["num_global_support_dofs"] = basix_element.num_global_support_dofs
    ir["num_element_support_dofs"] = basix_element.dim - ir["num_global_support_dofs"]

    return ir_dofmap(**ir)


def _compute_integral_ir(form_data, form_index, element_numbers, integral_names,
                         finite_element_names, parameters, visualise):
    """Compute intermediate represention for form integrals."""
    _entity_types = {
        "cell": "cell",
        "exterior_facet": "facet",
        "interior_facet": "facet",
        "vertex": "vertex",
        "custom": "cell"
    }

    # Iterate over groups of integrals
    irs = []
    for itg_data_index, itg_data in enumerate(form_data.integral_data):

        logger.info(f"Computing IR for integral in integral group {itg_data_index}")

        # Compute representation
        entitytype = _entity_types[itg_data.integral_type]
        cell = itg_data.domain.ufl_cell()
        cellname = cell.cellname()
        tdim = cell.topological_dimension()
        assert all(tdim == itg.ufl_domain().topological_dimension() for itg in itg_data.integrals)

        ir = {
            "integral_type": itg_data.integral_type,
            "subdomain_id": itg_data.subdomain_id,
            "rank": form_data.rank,
            "geometric_dimension": form_data.geometric_dimension,
            "topological_dimension": tdim,
            "entitytype": entitytype,
            "num_facets": cell.num_facets(),
            "num_vertices": cell.num_vertices(),
            "enabled_coefficients": itg_data.enabled_coefficients,
            "cell_shape": cellname,
            "coordinate_element": finite_element_names[itg_data.domain.ufl_coordinate_element()]
        }

        # Get element space dimensions
        unique_elements = element_numbers.keys()
        ir["element_dimensions"] = {
            ufl_element: create_element(ufl_element).dim
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

            if scheme == "custom":
                points = md["quadrature_points"]
                weights = md["quadrature_weights"]
            elif scheme == "vertex":
                # FIXME: Could this come from basix?

                # The vertex scheme, i.e., averaging the function value in the
                # vertices and multiplying with the simplex volume, is only of
                # order 1 and inferior to other generic schemes in terms of
                # error reduction. Equation systems generated with the vertex
                # scheme have some properties that other schemes lack, e.g., the
                # mass matrix is a simple diagonal matrix. This may be
                # prescribed in certain cases.
                if degree > 1:
                    warnings.warn(
                        "Explicitly selected vertex quadrature (degree 1), but requested degree is {}.".
                        format(degree))
                if cellname == "tetrahedron":
                    points, weights = (numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
                                                    [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
                                       numpy.array([1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0]))
                elif cellname == "triangle":
                    points, weights = (numpy.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
                                       numpy.array([1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0]))
                elif cellname == "interval":
                    # Trapezoidal rule
                    return (numpy.array([[0.0], [1.0]]), numpy.array([1.0 / 2.0, 1.0 / 2.0]))
            else:
                points, weights = create_quadrature_points_and_weights(
                    integral_type, cell, degree, scheme)

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
            _offset += numpy.product(constant.ufl_shape, dtype=int)

        ir["original_constant_offsets"] = original_constant_offsets

        ir["precision"] = itg_data.metadata["precision"]

        # Create map from number of quadrature points -> integrand
        integrands = {rule: integral.integrand() for rule, integral in sorted_integrals.items()}

        # Build more specific intermediate representation
        integral_ir = compute_integral_ir(itg_data.domain.ufl_cell(), itg_data.integral_type,
                                          ir["entitytype"], integrands, ir["tensor_shape"],
                                          parameters, visualise)

        ir.update(integral_ir)

        # Fetch name
        ir["name"] = integral_names[(form_index, itg_data_index)]

        irs.append(ir_integral(**ir))

    return irs


def _compute_form_ir(form_data, form_id, prefix, form_names, integral_names, element_numbers, finite_element_names,
                     dofmap_names, object_names):
    """Compute intermediate representation of form."""
    logger.info(f"Computing IR for form {form_id}")

    # Store id
    ir = {"id": form_id}

    # Compute common data
    ir["name"] = form_names[form_id]

    ir["signature"] = form_data.original_form.signature()

    ir["rank"] = len(form_data.original_form.arguments())
    ir["num_coefficients"] = len(form_data.reduced_coefficients)
    ir["num_constants"] = len(form_data.original_form.constants())

    ir["coefficient_names"] = [object_names.get(id(obj), f"w{j}")
                               for j, obj in enumerate(form_data.reduced_coefficients)]

    ir["constant_names"] = [object_names.get(id(obj), f"c{j}")
                            for j, obj in enumerate(form_data.original_form.constants())]

    ir["original_coefficient_position"] = form_data.original_coefficient_positions

    ir["finite_elements"] = [
        finite_element_names[e]
        for e in form_data.argument_elements + form_data.coefficient_elements
    ]

    ir["dofmaps"] = [
        dofmap_names[e] for e in form_data.argument_elements + form_data.coefficient_elements
    ]

    fs = {}
    for function in form_data.original_form.arguments() + tuple(form_data.reduced_coefficients):
        name = object_names.get(id(function), str(function))
        el = function.ufl_function_space().ufl_element()
        cmap = function.ufl_function_space().ufl_domain().ufl_coordinate_element()
        # Default point spacing for CoordinateElement is equispaced
        if cmap.variant() is None:
            cmap._sub_element._variant = "equispaced"
        basix_cmap = create_element(cmap)
        family = cmap.family()
        degree = cmap.degree()
        fs[name] = (finite_element_names[el], dofmap_names[el], family, degree,
                    basix_cmap.cell_type, basix_cmap.lagrange_variant)

    form_name = object_names.get(id(form_data.original_form), form_id)

    ir["function_spaces"] = fs
    ir["name_from_uflfile"] = f"form_{prefix}_{form_name}"

    # Store names of integrals and subdomain_ids for this form, grouped by integral types
    # Since form points to all integrals it contains, it has to know their names
    # for codegen phase
    ir["integral_names"] = {}
    ir["subdomain_ids"] = {}
    ufcx_integral_types = ("cell", "exterior_facet", "interior_facet")
    for integral_type in ufcx_integral_types:
        ir["subdomain_ids"][integral_type] = []
        ir["integral_names"][integral_type] = []

        for itg_index, itg_data in enumerate(form_data.integral_data):
            if (itg_data.integral_type == integral_type):
                if itg_data.subdomain_id == "otherwise":
                    # UFL is using "otherwise" for default integrals (over whole mesh)
                    # but FFCx needs integers, so otherwise = -1
                    if len(ir["subdomain_ids"][integral_type]) > 0 and ir["subdomain_ids"][integral_type][0] == -1:
                        raise ValueError("Only one default ('otherwise') integral allowed.")

                    # Put default integral as first
                    ir["subdomain_ids"][integral_type] = [-1] + ir["subdomain_ids"][integral_type]
                    ir["integral_names"][integral_type] = [
                        integral_names[(form_id, itg_index)]] + ir["integral_names"][integral_type]
                elif itg_data.subdomain_id < 0:
                    raise ValueError("Integral subdomain ID must be non-negative.")
                else:
                    assert isinstance(itg_data.subdomain_id, int)
                    ir["subdomain_ids"][integral_type] += [itg_data.subdomain_id]
                    ir["integral_names"][integral_type] += [integral_names[(form_id, itg_index)]]

    return ir_form(**ir)


def _compute_expression_ir(expression, index, prefix, analysis, parameters, visualise, object_names,
                           finite_element_names, dofmap_names):
    """Compute intermediate representation of expression."""
    logger.info(f"Computing IR for expression {index}")

    # Compute representation
    ir = {}

    original_expression = (expression[2], expression[1])

    ir["name"] = naming.expression_name(original_expression, prefix)

    original_expression = expression[2]
    points = expression[1]
    expression = expression[0]

    try:
        cell = expression.ufl_domain().ufl_cell()
    except AttributeError:
        # This case corresponds to a spatially constant expression
        # without any dependencies
        cell = None

    # Prepare dimensions of all unique element in expression, including
    # elements for arguments, coefficients and coordinate mappings
    ir["element_dimensions"] = {
        ufl_element: create_element(ufl_element).dim
        for ufl_element in analysis.unique_elements
    }

    # Extract dimensions for elements of arguments only
    arguments = ufl.algorithms.extract_arguments(expression)
    argument_elements = tuple(f.ufl_function_space().ufl_element() for f in arguments)
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

    ir["coefficient_names"] = [object_names.get(id(obj), f"w{j}")
                               for j, obj in enumerate(coefficients)]

    ir["constant_names"] = [object_names.get(id(obj), f"c{j}")
                            for j, obj in enumerate(ufl.algorithms.analysis.extract_constants(expression))]

    fs = {}
    for function in tuple(original_coefficients) + tuple(arguments):
        name = object_names.get(id(function), str(function))
        el = function.ufl_function_space().ufl_element()
        cmap = function.ufl_function_space().ufl_domain().ufl_coordinate_element()
        family = cmap.family()
        degree = cmap.degree()
        fs[name] = (finite_element_names[el], dofmap_names[el], family, degree)

    expression_name = object_names.get(id(original_expression), index)

    ir["function_spaces"] = fs
    ir["name_from_uflfile"] = f"expression_{prefix}_{expression_name}"

    if len(argument_elements) > 1:
        raise RuntimeError("Expression with more than one Argument not implemented.")

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
        _offset += numpy.product(constant.ufl_shape, dtype=int)

    ir["original_constant_offsets"] = original_constant_offsets

    ir["points"] = points

    weights = numpy.array([1.0] * points.shape[0])
    rule = QuadratureRule(points, weights)
    integrands = {rule: expression}

    if cell is None:
        assert len(ir["original_coefficient_positions"]) == 0 and len(ir["original_constant_offsets"]) == 0

    expression_ir = compute_integral_ir(cell, ir["integral_type"], ir["entitytype"], integrands, tensor_shape,
                                        parameters, visualise)

    ir.update(expression_ir)

    return ir_expression(**ir)
