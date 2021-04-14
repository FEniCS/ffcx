# Copyright (C) 2009-2020 Anders Logg, Martin Sandve Alnæs, Marie E. Rognes,
# Kristian B. Oelgaard, Matthew W. Scroggs, Chris Richardson, and others
#
# This file is part of FFCX.(https://www.fenicsproject.org)
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
from ffcx.basix_interface import create_basix_element
from ffcx.ir.integral import compute_integral_ir
from ffcx.ir.representationutils import (QuadratureRule,
                                         create_quadrature_points_and_weights)
from ufl.classes import Integral
from ufl.sorting import sorted_expr_sum

logger = logging.getLogger("ffcx")

# List of supported integral types
ufc_integral_types = ("cell", "exterior_facet", "interior_facet", "vertex", "custom")

ir_form = namedtuple('ir_form', [
    'id', 'prefix', 'name', 'signature', 'rank', 'num_coefficients', 'num_constants',
    'name_from_uflfile', 'function_spaces', 'original_coefficient_position',
    'coefficient_names', 'constant_names', 'create_coordinate_mapping', 'create_finite_element',
    'create_dofmap', 'create_cell_integral', 'get_cell_integral_ids', 'create_exterior_facet_integral',
    'get_exterior_facet_integral_ids', 'create_interior_facet_integral',
    'get_interior_facet_integral_ids', 'create_vertex_integral', 'get_vertex_integral_ids',
    'create_custom_integral', 'get_custom_integral_ids'])
ir_element = namedtuple('ir_element', [
    'id', 'name', 'signature', 'cell_shape', 'topological_dimension',
    'geometric_dimension', 'space_dimension', 'value_shape', 'reference_value_shape', 'degree',
    'family', 'num_sub_elements', 'block_size', 'create_sub_element',
    'entity_dofs', 'base_transformations',
    'needs_transformation_data', 'interpolation_is_identity'])
ir_dofmap = namedtuple('ir_dofmap', [
    'id', 'name', 'signature', 'num_global_support_dofs', 'num_element_support_dofs', 'num_entity_dofs',
    'tabulate_entity_dofs', 'base_transformations', 'num_sub_dofmaps', 'create_sub_dofmap', 'block_size'])
ir_coordinate_map = namedtuple('ir_coordinate_map', [
    'id', 'prefix', 'name', 'signature', 'cell_shape', 'topological_dimension', 'geometric_dimension',
    'coordinate_element_degree', 'coordinate_element_family', 'scalar_dofmap_name'])
ir_integral = namedtuple('ir_integral', [
    'integral_type', 'subdomain_id', 'rank', 'geometric_dimension', 'topological_dimension', 'entitytype',
    'num_facets', 'num_vertices', 'enabled_coefficients', 'element_dimensions',
    'element_ids', 'tensor_shape', 'coefficient_numbering', 'coefficient_offsets',
    'original_constant_offsets', 'params', 'cell_shape', 'unique_tables', 'unique_table_types',
    'table_dofmaps', 'table_dof_base_transformations', 'integrand', 'name', 'precision',
    'table_needs_transformation_data', 'needs_transformation_data'])
ir_evaluate_dof = namedtuple('ir_evaluate_dof', [
    'mappings', 'reference_value_size', 'physical_value_size', 'geometric_dimension',
    'topological_dimension', 'dofs', 'cell_shape'])
ir_expression = namedtuple('ir_expression', [
    'name', 'element_dimensions', 'params', 'unique_tables', 'unique_table_types', 'integrand',
    'table_dofmaps', 'table_dof_base_transformations', 'coefficient_numbering', 'coefficient_offsets',
    'integral_type', 'entitytype', 'tensor_shape', 'expression_shape', 'original_constant_offsets',
    'original_coefficient_positions', 'points', 'table_needs_transformation_data', 'needs_transformation_data'])

ir_data = namedtuple('ir_data', ['elements', 'dofmaps', 'coordinate_mappings', 'integrals', 'forms', 'expressions'])


def compute_ir(analysis: namedtuple, object_names, prefix, parameters, visualise):
    """Compute intermediate representation.

    """

    logger.info(79 * "*")
    logger.info("Compiler stage 2: Computing intermediate representation of objects")
    logger.info(79 * "*")

    # Compute object names
    # NOTE: This is done here for performance reasons, because repeated calls
    # within each IR computation would be expensive due to UFL signature computations
    finite_element_names = {e: naming.finite_element_name(e, prefix) for e in analysis.unique_elements}
    dofmap_names = {e: naming.dofmap_name(e, prefix) for e in analysis.unique_elements}
    coordinate_mapping_names = {cmap: naming.coordinate_map_name(
        cmap, prefix) for cmap in analysis.unique_coordinate_elements}
    integral_names = {}
    for fd_index, fd in enumerate(analysis.form_data):
        for itg_index, itg_data in enumerate(fd.integral_data):
            integral_names[(fd_index, itg_index)] = naming.integral_name(itg_data.integral_type, fd.original_form,
                                                                         fd_index, itg_data.subdomain_id)

    ir_elements = [
        _compute_element_ir(e, analysis.element_numbers, finite_element_names, parameters["epsilon"])
        for e in analysis.unique_elements
    ]

    ir_dofmaps = [
        _compute_dofmap_ir(e, analysis.element_numbers, dofmap_names) for e in analysis.unique_elements
    ]

    ir_coordinate_mappings = [
        _compute_coordinate_mapping_ir(e, prefix, analysis.element_numbers,
                                       coordinate_mapping_names, dofmap_names, finite_element_names)
        for e in analysis.unique_coordinate_elements
    ]

    irs = [
        _compute_integral_ir(fd, i, prefix, analysis.element_numbers, integral_names, parameters, visualise)
        for (i, fd) in enumerate(analysis.form_data)
    ]
    ir_integrals = list(itertools.chain(*irs))

    ir_forms = [
        _compute_form_ir(fd, i, prefix, analysis.element_numbers, finite_element_names,
                         dofmap_names, coordinate_mapping_names, object_names)
        for (i, fd) in enumerate(analysis.form_data)
    ]

    ir_expressions = [_compute_expression_ir(expr, i, prefix, analysis, parameters, visualise)
                      for i, expr in enumerate(analysis.expressions)]

    return ir_data(elements=ir_elements, dofmaps=ir_dofmaps,
                   coordinate_mappings=ir_coordinate_mappings,
                   integrals=ir_integrals, forms=ir_forms,
                   expressions=ir_expressions)


def _compute_element_ir(ufl_element, element_numbers, finite_element_names, epsilon):
    """Compute intermediate representation of element."""

    logger.info(f"Computing IR for element {ufl_element}")

    # Create basix elements
    basix_element = create_basix_element(ufl_element)
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
    ir["degree"] = ufl_element.degree()
    ir["family"] = basix_element.family_name
    ir["value_shape"] = ufl_element.value_shape()
    ir["reference_value_shape"] = ufl_element.reference_value_shape()

    ir["num_sub_elements"] = ufl_element.num_sub_elements()
    ir["create_sub_element"] = [finite_element_names[e] for e in ufl_element.sub_elements()]

    if hasattr(basix_element, "block_size"):
        ir["block_size"] = basix_element.block_size
        ufl_element = ufl_element.sub_elements()[0]
        basix_element = create_basix_element(ufl_element)
    else:
        ir["block_size"] = 1

    im = basix_element.interpolation_matrix
    if im.shape[0] == im.shape[1] and numpy.allclose(im, numpy.identity(im.shape[0])):
        ir["interpolation_is_identity"] = 1
    else:
        ir["interpolation_is_identity"] = 0

    ir["base_transformations"] = basix_element.base_transformations
    ir["needs_transformation_data"] = 0
    for p in basix_element.base_transformations:
        if not numpy.allclose(p, numpy.identity(len(p))):
            ir["needs_transformation_data"] = 1

    ir["entity_dofs"] = basix_element.entity_dof_numbers

    return ir_element(**ir)


def _compute_dofmap_ir(ufl_element, element_numbers, dofmap_names):
    """Compute intermediate representation of dofmap."""

    logger.info(f"Computing IR for dofmap of {ufl_element}")

    # Create basix elements
    basix_element = create_basix_element(ufl_element)

    # Store id
    ir = {"id": element_numbers[ufl_element]}
    ir["name"] = dofmap_names[ufl_element]

    # Compute data for each function
    ir["signature"] = "FFCX dofmap for " + repr(ufl_element)
    ir["create_sub_dofmap"] = [dofmap_names[e] for e in ufl_element.sub_elements()]
    ir["num_sub_dofmaps"] = ufl_element.num_sub_elements()

    if hasattr(basix_element, "block_size"):
        ir["block_size"] = basix_element.block_size
        basix_element = basix_element.sub_element
    else:
        ir["block_size"] = 1

    ir["base_transformations"] = basix_element.base_transformations

    # Precompute repeatedly used items
    for i in basix_element.entity_dofs:
        if max(i) != min(i):
            raise RuntimeError("Elements with different numbers of DOFs on subentities of the same dimension"
                               " are not yet supported in FFCx.")
    num_dofs_per_entity = [i[0] for i in basix_element.entity_dofs]

    ir["num_entity_dofs"] = num_dofs_per_entity
    ir["tabulate_entity_dofs"] = (basix_element.entity_dof_numbers, num_dofs_per_entity)

    ir["num_global_support_dofs"] = basix_element.num_global_support_dofs
    ir["num_element_support_dofs"] = basix_element.dim - ir["num_global_support_dofs"]

    return ir_dofmap(**ir)


_midpoints = {
    "interval": (0.5, ),
    "triangle": (1.0 / 3.0, 1.0 / 3.0),
    "tetrahedron": (0.25, 0.25, 0.25),
    "quadrilateral": (0.5, 0.5),
    "hexahedron": (0.5, 0.5, 0.5),
}


def cell_midpoint(cell):
    # TODO: Is this defined somewhere more central where we can get it from?
    return _midpoints[cell.cellname()]


def _tabulate_coordinate_mapping_basis(ufl_element):
    # TODO: Move this function to a table generation module?

    # Get scalar element, assuming coordinates are represented
    # with a VectorElement of scalar subelements
    selement = ufl_element.sub_elements()[0]

    basix_element = create_basix_element(selement)
    cell = selement.cell()
    tdim = cell.topological_dimension()

    tables = {}

    # Get points
    origin = (0.0, ) * tdim
    midpoint = cell_midpoint(cell)

    # Tabulate basis
    t0 = basix_element.tabulate(1, [origin])
    tm = basix_element.tabulate(1, [midpoint])

    # Get basis values at cell origin
    tables["x0"] = t0[0][:, 0]

    # Get basis values at cell midpoint
    tables["xm"] = tm[0][:, 0]

    # Get basis derivative values at cell origin
    tables["J0"] = numpy.asarray([t0[d][:, 0] for d in range(1, 1 + tdim)])

    # Get basis derivative values at cell midpoint
    tables["Jm"] = numpy.asarray([tm[d][:, 0] for d in range(1, 1 + tdim)])

    return tables


def _compute_coordinate_mapping_ir(ufl_coordinate_element,
                                   prefix,
                                   element_numbers,
                                   coordinate_mapping_names,
                                   dofmap_names,
                                   finite_element_names):
    """Compute intermediate representation of coordinate mapping."""

    logger.info(f"Computing IR for coordinate mapping {ufl_coordinate_element}")

    cell = ufl_coordinate_element.cell()
    cellname = cell.cellname()

    assert ufl_coordinate_element.value_shape() == (cell.geometric_dimension(), )

    # Compute element values
    tables = _tabulate_coordinate_mapping_basis(ufl_coordinate_element)

    # Store id
    ir = {"id": element_numbers[ufl_coordinate_element]}
    ir["prefix"] = prefix
    ir["name"] = coordinate_mapping_names[ufl_coordinate_element]

    # Compute data for each function
    ir["signature"] = "FFCX coordinate_mapping from " + repr(ufl_coordinate_element)
    ir["cell_shape"] = cellname
    ir["topological_dimension"] = cell.topological_dimension()
    ir["geometric_dimension"] = ufl_coordinate_element.value_size()

    # NB! The entries below breaks the pattern of using ir keywords == code keywords,
    # which I personally don't find very useful anyway (martinal).

    basix_element = create_basix_element(ufl_coordinate_element)

    # Store tables and other coordinate element data
    ir["coordinate_element_degree"] = ufl_coordinate_element.degree()
    ir["coordinate_element_family"] = basix_element.family_name

    # Get classnames for finite element and dofmap of scalar subelement
    scalar_element = ufl_coordinate_element.sub_elements()[0]
    ir["scalar_dofmap_name"] = dofmap_names[scalar_element]

    return ir_coordinate_map(**ir)


def _compute_integral_ir(form_data, form_index, prefix, element_numbers, integral_names,
                         parameters, visualise):
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
            "cell_shape": cellname
        }

        # Get element space dimensions
        unique_elements = element_numbers.keys()
        ir["element_dimensions"] = {
            ufl_element: create_basix_element(ufl_element).dim
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


def _compute_form_ir(form_data, form_id, prefix, element_numbers, finite_element_names,
                     dofmap_names, coordinate_mapping_names, object_names):
    """Compute intermediate representation of form."""

    logger.info(f"Computing IR for form {form_id}")

    # Store id
    ir = {"id": form_id}

    # Storing prefix here for reconstruction of classnames on code
    # generation side
    ir["prefix"] = prefix

    # Compute common data
    ir["name"] = naming.form_name(form_data.original_form, form_id)

    ir["signature"] = form_data.original_form.signature()

    ir["rank"] = len(form_data.original_form.arguments())
    ir["num_coefficients"] = len(form_data.reduced_coefficients)
    ir["num_constants"] = len(form_data.original_form.constants())

    ir["coefficient_names"] = [object_names.get(id(obj), "w%d" % j)
                               for j, obj in enumerate(form_data.reduced_coefficients)]

    ir["constant_names"] = [object_names.get(id(obj), "c%d" % j)
                            for j, obj in enumerate(form_data.original_form.constants())]

    ir["original_coefficient_position"] = form_data.original_coefficient_positions

    ir["create_coordinate_mapping"] = [
        coordinate_mapping_names[e] for e in form_data.coordinate_elements
    ]
    ir["create_finite_element"] = [
        finite_element_names[e]
        for e in form_data.argument_elements + form_data.coefficient_elements
    ]
    ir["create_dofmap"] = [
        dofmap_names[e] for e in form_data.argument_elements + form_data.coefficient_elements
    ]

    fs = {}
    for function in form_data.original_form.arguments() + tuple(form_data.reduced_coefficients):
        name = object_names.get(id(function), str(function))
        el = function.ufl_element()
        cmap = function.ufl_function_space().ufl_domain().ufl_coordinate_element()
        fs[name] = (finite_element_names[el], dofmap_names[el], coordinate_mapping_names[cmap])

    form_name = object_names.get(id(form_data.original_form), form_id)

    ir["function_spaces"] = fs
    ir["name_from_uflfile"] = f"form_{prefix}_{form_name}"

    # Create integral ids and names using form prefix (integrals are
    # always generated as part of form so don't get their own prefix)
    for integral_type in ufc_integral_types:
        irdata = _create_foo_integral(prefix, form_id, integral_type, form_data)
        ir[f"create_{integral_type}_integral"] = irdata
        ir[f"get_{integral_type}_integral_ids"] = irdata

    return ir_form(**ir)


def _compute_expression_ir(expression, index, prefix, analysis, parameters, visualise):

    logger.info(f"Computing IR for expression {index}")

    # Compute representation
    ir = {}

    original_expression = (expression[2], expression[1])
    sig = naming.compute_signature([original_expression], "", parameters)
    ir["name"] = "expression_{!s}".format(sig)

    original_expression = expression[2]
    points = expression[1]
    expression = expression[0]

    try:
        cell = expression.ufl_domain().ufl_cell()
    except AttributeError:
        # This case corresponds to a spatially constant expression without any dependencies
        cell = None

    # Prepare dimensions of all unique element in expression, including
    # elements for arguments, coefficients and coordinate mappings
    ir["element_dimensions"] = {
        ufl_element: create_basix_element(ufl_element).dim
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


def _create_foo_integral(prefix, form_id, integral_type, form_data):
    """Compute intermediate representation of create_foo_integral."""
    subdomain_ids = []
    classnames = []
    itg_data = [itg_data for itg_data in form_data.integral_data
                if (itg_data.integral_type == integral_type and itg_data.subdomain_id == "otherwise")]

    if len(itg_data) > 1:
        raise RuntimeError("Expecting at most one default integral of each type.")
    elif len(itg_data) == 1:
        subdomain_ids += [-1]
        classnames += [naming.integral_name(integral_type, form_data.original_form,
                                            form_id, "otherwise")]

    for itg_data in form_data.integral_data:
        if isinstance(itg_data.subdomain_id, int):
            if itg_data.subdomain_id < 0:
                raise ValueError(f"Integral subdomain ID must be non-negative, not {itg_data.subdomain_id}")
            if (itg_data.integral_type == integral_type):
                subdomain_ids += [itg_data.subdomain_id]
                classnames += [naming.integral_name(integral_type, form_data.original_form,
                                                    form_id, itg_data.subdomain_id)]

    return subdomain_ids, classnames
