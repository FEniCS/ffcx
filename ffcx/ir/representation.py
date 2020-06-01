# Copyright (C) 2009-2017 Anders Logg, Martin Sandve Alnæs, Marie E. Rognes,
# Kristian B. Oelgaard, and others
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
from collections import namedtuple

import numpy

import ufl
from ffcx import naming
from ffcx.fiatinterface import (EnrichedElement, FlattenedDimensions,
                                MixedElement, QuadratureElement, SpaceOfReals,
                                create_element)
from ffcx.ir import dof_permutations
from ffcx.ir.representationutils import (QuadratureRule,
                                         create_quadrature_points_and_weights)
from ffcx.ir.integral import compute_integral_ir
from ufl.sorting import sorted_expr_sum
from ufl.classes import Integral
from FIAT.hdiv_trace import HDivTrace

logger = logging.getLogger(__name__)

# List of supported integral types
ufc_integral_types = ("cell", "exterior_facet", "interior_facet", "vertex", "custom")

ir_form = namedtuple('ir_form', ['id', 'prefix', 'name', 'signature', 'rank',
                                 'num_coefficients', 'num_constants', 'name_from_uflfile',
                                 'function_spaces',
                                 'original_coefficient_position',
                                 'coefficient_names', 'constant_names',
                                 'create_coordinate_mapping', 'create_finite_element',
                                 'create_dofmap', 'create_cell_integral',
                                 'get_cell_integral_ids', 'create_exterior_facet_integral',
                                 'get_exterior_facet_integral_ids', 'create_interior_facet_integral',
                                 'get_interior_facet_integral_ids', 'create_vertex_integral',
                                 'get_vertex_integral_ids', 'create_custom_integral',
                                 'get_custom_integral_ids'])
ir_element = namedtuple('ir_element', ['id', 'name', 'signature', 'cell_shape',
                                       'topological_dimension',
                                       'geometric_dimension', 'space_dimension', 'value_shape',
                                       'reference_value_shape', 'degree', 'family', 'evaluate_basis',
                                       'evaluate_dof', 'tabulate_dof_coordinates', 'num_sub_elements',
                                       'base_permutations', 'dof_reflection_entities',
                                       'create_sub_element', 'dof_types', 'entity_dofs'])
ir_dofmap = namedtuple('ir_dofmap', ['id', 'name', 'signature', 'num_global_support_dofs',
                                     'num_element_support_dofs', 'num_entity_dofs',
                                     'tabulate_entity_dofs', 'base_permutations', 'dof_reflection_entities',
                                     'num_sub_dofmaps', 'create_sub_dofmap', 'dof_types'])
ir_coordinate_map = namedtuple('ir_coordinate_map', ['id', 'prefix', 'name', 'signature', 'cell_shape',
                                                     'topological_dimension',
                                                     'geometric_dimension',
                                                     'compute_physical_coordinates',
                                                     'compute_reference_coordinates', 'compute_jacobians',
                                                     'compute_jacobian_determinants',
                                                     'compute_jacobian_inverses', 'compute_geometry', 'tables',
                                                     'coordinate_element_degree', 'num_scalar_coordinate_element_dofs',
                                                     'coordinate_finite_element_classname',
                                                     'scalar_coordinate_finite_element_classname',
                                                     'scalar_dofmap_name', 'is_affine'])
ir_integral = namedtuple('ir_integral', ['integral_type', 'subdomain_id',
                                         'rank', 'geometric_dimension', 'topological_dimension',
                                         'entitytype', 'num_facets', 'num_vertices', 'needs_oriented',
                                         'enabled_coefficients', 'element_dimensions', 'element_ids',
                                         'tensor_shape', 'coefficient_numbering',
                                         'coefficient_offsets', 'original_constant_offsets', 'params', 'cell_shape',
                                         'unique_tables', 'unique_table_types', 'table_dofmaps',
                                         'table_dof_face_tangents', 'table_dof_reflection_entities',
                                         'integrand', 'name', 'precision'])
ir_tabulate_dof_coordinates = namedtuple('ir_tabulate_dof_coordinates', ['tdim', 'gdim', 'points', 'cell_shape'])
ir_evaluate_dof = namedtuple('ir_evaluate_dof', ['mappings', 'reference_value_size', 'physical_value_size',
                                                 'geometric_dimension', 'topological_dimension', 'dofs',
                                                 'physical_offsets', 'cell_shape'])
ir_expression = namedtuple('ir_expression', ['name', 'element_dimensions', 'params', 'unique_tables',
                                             'unique_table_types', 'integrand', 'table_dofmaps',
                                             'table_dof_face_tangents', 'table_dof_reflection_entities',
                                             'coefficient_numbering', 'coefficient_offsets',
                                             'integral_type', 'entitytype', 'tensor_shape', 'expression_shape',
                                             'original_constant_offsets', 'original_coefficient_positions', 'points'])

ir_data = namedtuple('ir_data', ['elements', 'dofmaps', 'coordinate_mappings', 'integrals', 'forms', 'expressions'])


def compute_ir(analysis: namedtuple, object_names, prefix, parameters, visualise):
    """Compute intermediate representation.

    """
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

    logger.info("Computing IR for element {}".format(ufl_element))

    # Create FIAT element
    fiat_element = create_element(ufl_element)
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
    ir["space_dimension"] = fiat_element.space_dimension()
    ir["value_shape"] = ufl_element.value_shape()
    ir["reference_value_shape"] = ufl_element.reference_value_shape()

    ir["degree"] = ufl_element.degree()
    ir["family"] = ufl_element.family()

    ir["evaluate_basis"] = _evaluate_basis(ufl_element, fiat_element, epsilon)
    ir["evaluate_dof"] = _evaluate_dof(ufl_element, fiat_element)
    ir["tabulate_dof_coordinates"] = _tabulate_dof_coordinates(ufl_element, fiat_element)
    ir["num_sub_elements"] = ufl_element.num_sub_elements()
    ir["create_sub_element"] = [finite_element_names[e] for e in ufl_element.sub_elements()]

    ir["base_permutations"] = dof_permutations.base_permutations(ufl_element)
    ir["dof_reflection_entities"] = dof_permutations.reflection_entities(ufl_element)

    ir["dof_types"] = [i.functional_type for i in fiat_element.dual_basis()]
    ir["entity_dofs"] = fiat_element.entity_dofs()

    return ir_element(**ir)


def _compute_dofmap_ir(ufl_element, element_numbers, dofmap_names):
    """Compute intermediate representation of dofmap."""

    logger.info("Computing IR for dofmap of {}".format(ufl_element))

    # Create FIAT element
    fiat_element = create_element(ufl_element)

    # Precompute repeatedly used items
    num_dofs_per_entity = _num_dofs_per_entity(fiat_element)
    entity_dofs = fiat_element.entity_dofs()

    # Store id
    ir = {"id": element_numbers[ufl_element]}
    ir["name"] = dofmap_names[ufl_element]

    # Compute data for each function
    ir["signature"] = "FFCX dofmap for " + repr(ufl_element)
    ir["num_global_support_dofs"] = _num_global_support_dofs(fiat_element)
    ir["num_element_support_dofs"] = fiat_element.space_dimension() - ir["num_global_support_dofs"]
    ir["num_entity_dofs"] = num_dofs_per_entity
    ir["tabulate_entity_dofs"] = (entity_dofs, num_dofs_per_entity)
    ir["num_sub_dofmaps"] = ufl_element.num_sub_elements()
    ir["create_sub_dofmap"] = [dofmap_names[e] for e in ufl_element.sub_elements()]
    ir["dof_types"] = [i.functional_type for i in fiat_element.dual_basis()]
    ir["base_permutations"] = dof_permutations.base_permutations(ufl_element)
    ir["dof_reflection_entities"] = dof_permutations.reflection_entities(ufl_element)

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

    fiat_element = create_element(selement)
    cell = selement.cell()
    tdim = cell.topological_dimension()

    tables = {}

    # Get points
    origo = (0.0, ) * tdim
    midpoint = cell_midpoint(cell)

    # Tabulate basis
    t0 = fiat_element.tabulate(1, [origo])
    tm = fiat_element.tabulate(1, [midpoint])

    # Get basis values at cell origo
    tables["x0"] = t0[(0, ) * tdim][:, 0]

    # Get basis values at cell midpoint
    tables["xm"] = tm[(0, ) * tdim][:, 0]

    # Single direction derivatives, e.g. [(1,0), (0,1)] in 2d
    derivatives = [(0, ) * i + (1, ) + (0, ) * (tdim - 1 - i) for i in range(tdim)]

    # Get basis derivative values at cell origo
    tables["J0"] = numpy.asarray([t0[d][:, 0] for d in derivatives])

    # Get basis derivative values at cell midpoint
    tables["Jm"] = numpy.asarray([tm[d][:, 0] for d in derivatives])

    return tables


def _compute_coordinate_mapping_ir(ufl_coordinate_element,
                                   prefix,
                                   element_numbers,
                                   coordinate_mapping_names,
                                   dofmap_names,
                                   finite_element_names):
    """Compute intermediate representation of coordinate mapping."""

    logger.info("Computing IR for coordinate mapping {}".format(ufl_coordinate_element))

    cell = ufl_coordinate_element.cell()
    cellname = cell.cellname()

    assert ufl_coordinate_element.value_shape() == (cell.geometric_dimension(), )

    # Compute element values via fiat element
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

    ir["compute_physical_coordinates"] = None  # currently unused, corresponds to function name
    ir["compute_reference_coordinates"] = None  # currently unused, corresponds to function name
    ir["compute_jacobians"] = None  # currently unused, corresponds to function name
    ir["compute_jacobian_determinants"] = None  # currently unused, corresponds to function name
    ir["compute_jacobian_inverses"] = None  # currently unused, corresponds to function name
    ir["compute_geometry"] = None  # currently unused, corresponds to function name

    # NB! The entries below breaks the pattern of using ir keywords == code keywords,
    # which I personally don't find very useful anyway (martinal).

    # Store tables and other coordinate element data
    ir["tables"] = tables
    ir["coordinate_element_degree"] = ufl_coordinate_element.degree()
    ir["num_scalar_coordinate_element_dofs"] = tables["x0"].shape[0]
    ir["is_affine"] = ir["coordinate_element_degree"] == 1 and cellname in ("interval", "triangle", "tetrahedron")

    # Get classnames for coordinate element
    ir["coordinate_finite_element_classname"] = finite_element_names[ufl_coordinate_element]

    # Get classnames for finite element and dofmap of scalar subelement
    scalar_element = ufl_coordinate_element.sub_elements()[0]
    ir["scalar_coordinate_finite_element_classname"] = finite_element_names[scalar_element]
    ir["scalar_dofmap_name"] = dofmap_names[scalar_element]

    return ir_coordinate_map(**ir)


def _num_global_support_dofs(fiat_element):
    """Compute number of global support dofs."""
    if not isinstance(fiat_element, MixedElement):
        if isinstance(fiat_element, SpaceOfReals):
            return 1
        return 0
    num_reals = 0
    for e in fiat_element.elements():
        if isinstance(e, SpaceOfReals):
            num_reals += 1
    return num_reals


def _compute_integral_ir(form_data, form_index, prefix, element_numbers, integral_names,
                         parameters, visualise):
    """Compute intermediate represention for form integrals."""

    logger.info("Computing IR for integrals {}".format(form_data))

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
            "needs_oriented": form_needs_oriented_jacobian(form_data),
            "enabled_coefficients": itg_data.enabled_coefficients,
            "cell_shape": cellname
        }

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

    logger.info("Computing IR for form {}".format(form_data))

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
    ir["name_from_uflfile"] = "form_{}_{}".format(prefix, form_name)

    # Create integral ids and names using form prefix (integrals are
    # always generated as part of form so don't get their own prefix)
    for integral_type in ufc_integral_types:
        irdata = _create_foo_integral(prefix, form_id, integral_type, form_data)
        ir["create_{}_integral".format(integral_type)] = irdata
        ir["get_{}_integral_ids".format(integral_type)] = irdata

    return ir_form(**ir)


def _compute_expression_ir(expression, index, prefix, analysis, parameters, visualise):

    logger.info("Computing IR for expression {}".format(str(expression)))

    # Compute representation
    ir = {}

    original_expression = (expression[2], expression[1])
    sig = naming.compute_signature([original_expression], "", parameters)
    ir["name"] = "expression_{!s}".format(sig)

    original_expression = expression[2]
    points = expression[1]
    expression = expression[0]

    cell = expression.ufl_domain().ufl_cell()

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

    expression_ir = compute_integral_ir(cell, ir["integral_type"], ir["entitytype"], integrands, tensor_shape,
                                        parameters, visualise)

    ir.update(expression_ir)

    return ir_expression(**ir)


def _generate_reference_offsets(fiat_element, offset=0):
    """Generate offsets.

    I.e., value offset for each basis function relative to a reference
    element representation.

    """
    if isinstance(fiat_element, MixedElement):
        offsets = []
        for e in fiat_element.elements():
            offsets += _generate_reference_offsets(e, offset)
            # NB! This is the fiat element and therefore value_shape
            # means reference_value_shape
            offset += ufl.utils.sequences.product(e.value_shape())
        return offsets
    elif isinstance(fiat_element, EnrichedElement):
        offsets = []
        for e in fiat_element.elements():
            offsets += _generate_reference_offsets(e, offset)
        return offsets
    else:
        return [offset] * fiat_element.space_dimension()


def _generate_physical_offsets(ufl_element, offset=0):
    """Generate offsets.

    I.e., value offset for each basis function relative to a physical
    element representation.

    """
    cell = ufl_element.cell()
    gdim = cell.geometric_dimension()
    tdim = cell.topological_dimension()

    # Refer to reference if gdim == tdim. This is a hack to support more
    # stuff (in particular restricted elements)
    if gdim == tdim:
        return _generate_reference_offsets(create_element(ufl_element))

    if isinstance(ufl_element, ufl.MixedElement):
        offsets = []
        for e in ufl_element.sub_elements():
            offsets += _generate_physical_offsets(e, offset)
            # e is a ufl element, so value_size means the physical value size
            offset += e.value_size()
        return offsets
    elif isinstance(ufl_element, ufl.EnrichedElement):
        offsets = []
        for e in ufl_element._elements:  # TODO: Avoid private member access
            offsets += _generate_physical_offsets(e, offset)
        return offsets
    elif isinstance(ufl_element, ufl.FiniteElement):
        fiat_element = create_element(ufl_element)
        return [offset] * fiat_element.space_dimension()
    else:
        raise NotImplementedError("This element combination is not implemented")


def _generate_offsets(ufl_element, reference_offset=0, physical_offset=0):
    """Generate offsets.

    I.e., value offset for each basis function relative to a physical
    element representation.

    """
    if isinstance(ufl_element, ufl.MixedElement):
        offsets = []
        for e in ufl_element.sub_elements():
            offsets += _generate_offsets(e, reference_offset, physical_offset)
            # e is a ufl element, so value_size means the physical value size
            reference_offset += e.reference_value_size()
            physical_offset += e.value_size()
        return offsets
    elif isinstance(ufl_element, ufl.EnrichedElement):
        offsets = []
        for e in ufl_element._elements:  # TODO: Avoid private member access
            offsets += _generate_offsets(e, reference_offset, physical_offset)
        return offsets
    elif isinstance(ufl_element, ufl.FiniteElement):
        fiat_element = create_element(ufl_element)
        return [(reference_offset, physical_offset)] * fiat_element.space_dimension()
    else:
        # TODO: Support RestrictedElement, QuadratureElement,
        #       TensorProductElement, etc.!  and replace
        #       _generate_{physical|reference}_offsets with this
        #       function.
        raise NotImplementedError("This element combination is not implemented")


def _evaluate_dof(ufl_element, fiat_element):
    """Compute intermediate representation of evaluate_dof."""
    cell = ufl_element.cell()
    if fiat_element.is_nodal():
        dofs = [L.pt_dict for L in fiat_element.dual_basis()]
    else:
        dofs = [None] * fiat_element.space_dimension()

    return ir_evaluate_dof(mappings=fiat_element.mapping(),
                           reference_value_size=ufl_element.reference_value_size(),
                           physical_value_size=ufl_element.value_size(),
                           geometric_dimension=cell.geometric_dimension(),
                           topological_dimension=cell.topological_dimension(),
                           dofs=dofs,
                           physical_offsets=_generate_physical_offsets(ufl_element),
                           cell_shape=cell.cellname())


def _extract_elements(fiat_element):
    new_elements = []
    if isinstance(fiat_element, (MixedElement, EnrichedElement)):
        for e in fiat_element.elements():
            new_elements += _extract_elements(e)
    else:
        new_elements.append(fiat_element)
    return new_elements


def _evaluate_basis(ufl_element, fiat_element, epsilon):
    """Compute intermediate representation for evaluate_basis."""
    cell = ufl_element.cell()
    cellname = cell.cellname()

    # Handle Mixed and EnrichedElements by extracting 'sub' elements.
    elements = _extract_elements(fiat_element)
    physical_offsets = _generate_physical_offsets(ufl_element)
    reference_offsets = _generate_reference_offsets(fiat_element)
    mappings = fiat_element.mapping()

    # This function is evidently not implemented for TensorElements
    for e in elements:
        if (len(e.value_shape()) > 1) and (e.num_sub_elements() != 1):
            return "Function not supported/implemented for TensorElements."

    # Handle QuadratureElement, not supported because the basis is only
    # defined at the dof coordinates where the value is 1, so not very
    # interesting.
    for e in elements:
        if isinstance(e, QuadratureElement):
            return "Function not supported/implemented for QuadratureElement."
        if isinstance(e, FlattenedDimensions) and isinstance(e.element, QuadratureElement):
            # Case for quad/hex cell
            return "Function not supported/implemented for QuadratureElement."
        if isinstance(e, HDivTrace):
            return "Function not supported for Trace elements"

    # Initialise data with 'global' values.
    data = {
        "reference_value_size": ufl_element.reference_value_size(),
        "physical_value_size": ufl_element.value_size(),
        "cellname": cellname,
        "topological_dimension": cell.topological_dimension(),
        "geometric_dimension": cell.geometric_dimension(),
        "space_dimension": fiat_element.space_dimension(),
        "needs_oriented": element_needs_oriented_jacobian(fiat_element),
        "max_degree": max([e.degree() for e in elements])
    }

    # Loop element and space dimensions to generate dof data.
    dof = 0
    dofs_data = []
    for e in elements:
        num_components = ufl.utils.sequences.product(e.value_shape())
        if isinstance(e, FlattenedDimensions):
            # Tensor product element
            A = e.element.A
            B = e.element.B
            # Attach suitable coefficients to element
            if isinstance(A, FlattenedDimensions):
                # This is for hexahedral element
                ac = A.element.A.get_coeffs()
                bc = A.element.B.get_coeffs()
                ac = numpy.block([[w * ac for w in v] for v in bc])
                ad = A.element.A.dmats()
                bd = A.element.B.dmats()
                ai = numpy.eye(ad[0].shape[0])
                bi = numpy.eye(bd[0].shape[0])

                if len(bd) != 1:
                    raise NotImplementedError("Cannot create dmats")

                dmats = []
                for mat in ad:
                    dmats += [numpy.block([[w * mat for w in v] for v in bi])]
                dmats += [numpy.block([[w * ai for w in v] for v in bd[0]])]
                ad = dmats
            else:
                ac = A.get_coeffs()
                ad = A.dmats()
            bc = B.get_coeffs()
            bd = B.dmats()
            coeffs = numpy.block([[w * ac for w in v] for v in bc])
            num_expansion_members = coeffs.shape[0]
            ai = numpy.eye(ad[0].shape[0])
            bi = numpy.eye(bd[0].shape[0])

            if len(bd) != 1:
                raise NotImplementedError("Cannot create dmats")

            dmats = []
            for mat in ad:
                dmats += [numpy.block([[w * mat for w in v] for v in bi])]
            dmats += [numpy.block([[w * ai for w in v] for v in bd[0]])]

        else:
            coeffs = e.get_coeffs()
            dmats = e.dmats()
            num_expansion_members = e.get_num_members(e.degree())

        # Clamp dmats zeros
        dmats = numpy.asarray(dmats)
        dmats[numpy.where(numpy.isclose(dmats, 0.0, rtol=epsilon, atol=epsilon))] = 0.0

        # Extracted parts of dd below that are common for the element
        # here.  These dict entries are added to each dof_data dict for
        # each dof, because that's what the code generation
        # implementation expects.  If the code generation needs this
        # structure to be optimized in the future, we can store this
        # data for each subelement instead of for each dof.
        subelement_data = {
            "embedded_degree": e.degree(),
            "num_components": num_components,
            "dmats": dmats,
            "num_expansion_members": num_expansion_members,
        }
        value_rank = len(e.value_shape())

        for i in range(e.space_dimension()):
            if num_components == 1:
                coefficients = [coeffs[i]]
            elif value_rank == 1:
                # Handle coefficients for vector valued basis elements
                # [Raviart-Thomas, Brezzi-Douglas-Marini (BDM)].
                coefficients = [coeffs[i][c] for c in range(num_components)]
            elif value_rank == 2:
                # Handle coefficients for tensor valued basis elements.
                # [Regge]
                coefficients = [
                    coeffs[i][p][q] for p in range(e.value_shape()[0])
                    for q in range(e.value_shape()[1])
                ]
            else:
                raise RuntimeError("Unknown situation with num_components > 1")

            # Clamp coefficient zeros
            coefficients = numpy.asarray(coefficients)
            coefficients[numpy.where(numpy.isclose(coefficients, 0.0, rtol=epsilon,
                                                   atol=epsilon))] = 0.0

            dof_data = {
                "coeffs": coefficients,
                "mapping": mappings[dof],
                "physical_offset": physical_offsets[dof],
                "reference_offset": reference_offsets[dof],
            }
            # Still storing element data in dd to avoid rewriting dependent code
            dof_data.update(subelement_data)

            # This list will hold one dd dict for each dof
            dofs_data.append(dof_data)
            dof += 1

    data["dofs_data"] = dofs_data

    return data


def _tabulate_dof_coordinates(ufl_element, element):
    """Compute intermediate representation of tabulate_dof_coordinates."""
    if uses_integral_moments(element):
        return {}

    # Bail out if any dual basis member is missing (element is not
    # nodal), this is strictly not necessary but simpler
    if any(L is None for L in element.dual_basis()):
        return {}

    cell = ufl_element.cell()
    return ir_tabulate_dof_coordinates(
        tdim=cell.topological_dimension(),
        gdim=cell.geometric_dimension(),
        points=[sorted(L.pt_dict.keys())[0] for L in element.dual_basis()],
        cell_shape=cell.cellname())


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
                raise ValueError("Integral subdomain ID must be non-negative, not {}".format(itg_data.subdomain_id))
            if (itg_data.integral_type == integral_type):
                subdomain_ids += [itg_data.subdomain_id]
                classnames += [naming.integral_name(integral_type, form_data.original_form,
                                                    form_id, itg_data.subdomain_id)]

    return subdomain_ids, classnames


def all_elements(fiat_element):
    if isinstance(fiat_element, MixedElement):
        return fiat_element.elements()
    return [fiat_element]


def _num_dofs_per_entity(fiat_element):
    """Compute list of the number of dofs associated with a single mesh entity.

    Example: Lagrange of degree 3 on triangle: [1, 2, 1]

    """
    entity_dofs = fiat_element.entity_dofs()
    return [len(entity_dofs[e][0]) for e in sorted(entity_dofs.keys())]


def uses_integral_moments(fiat_element):
    """True if element uses integral moments for its degrees of freedom."""
    integrals = set(["IntegralMoment", "FrobeniusIntegralMoment"])
    tags = set([L.get_type_tag() for L in fiat_element.dual_basis() if L])
    return len(integrals & tags) > 0


def element_needs_oriented_jacobian(fiat_element):
    # Check whether this element needs an oriented jacobian (only
    # contravariant piolas seem to need it)
    return "contravariant piola" in fiat_element.mapping()


def form_needs_oriented_jacobian(form_data):
    # Check whether this form needs an oriented jacobian (only forms
    # involving contravariant piola mappings seem to need it)
    for ufl_element in form_data.unique_elements:
        element = create_element(ufl_element)
        if "contravariant piola" in element.mapping():
            return True
    return False
