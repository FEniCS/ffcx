# -*- coding: utf-8 -*-
# Copyright (C) 2009-2017 Anders Logg, Martin Sandve Alnæs, Marie E. Rognes,
# Kristian B. Oelgaard, and others
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Compiler stage 2: Code representation

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
from ffc import classname
from ffc.fiatinterface import (EnrichedElement, FlattenedDimensions,
                               MixedElement, QuadratureElement, SpaceOfReals,
                               create_element)
from FIAT.hdiv_trace import HDivTrace

logger = logging.getLogger(__name__)

# List of supported integral types
ufc_integral_types = ("cell", "exterior_facet", "interior_facet", "vertex", "custom")

ir_form = namedtuple('ir_form', ['id', 'prefix', 'classname', 'signature', 'rank',
                                 'num_coefficients', 'original_coefficient_position',
                                 'coefficient_names',
                                 'create_coordinate_finite_element', 'create_coordinate_dofmap',
                                 'create_coordinate_mapping', 'create_finite_element',
                                 'create_dofmap', 'create_cell_integral',
                                 'get_cell_integral_ids', 'create_exterior_facet_integral',
                                 'get_exterior_facet_integral_ids', 'create_interior_facet_integral',
                                 'get_interior_facet_integral_ids', 'create_vertex_integral',
                                 'get_vertex_integral_ids', 'create_custom_integral',
                                 'get_custom_integral_ids'])
ir_element = namedtuple('ir_element', ['id', 'classname', 'signature', 'cell_shape',
                                       'topological_dimension',
                                       'geometric_dimension', 'space_dimension', 'value_shape',
                                       'reference_value_shape', 'degree', 'family', 'evaluate_basis',
                                       'evaluate_dof', 'tabulate_dof_coordinates', 'num_sub_elements',
                                       'create_sub_element'])
ir_dofmap = namedtuple('ir_dofmap', ['id', 'classname', 'signature', 'num_global_support_dofs',
                                     'num_element_support_dofs', 'num_entity_dofs',
                                     'tabulate_entity_dofs',
                                     'num_sub_dofmaps', 'create_sub_dofmap'])
ir_coordinate_map = namedtuple('ir_coordinate_map', ['id', 'classname', 'signature', 'cell_shape',
                                                     'topological_dimension',
                                                     'geometric_dimension', 'create_coordinate_finite_element',
                                                     'create_coordinate_dofmap', 'compute_physical_coordinates',
                                                     'compute_reference_coordinates', 'compute_jacobians',
                                                     'compute_jacobian_determinants',
                                                     'compute_jacobian_inverses', 'compute_geometry', 'tables',
                                                     'coordinate_element_degree', 'num_scalar_coordinate_element_dofs',
                                                     'coordinate_finite_element_classname',
                                                     'scalar_coordinate_finite_element_classname'])
ir_integral = namedtuple('ir_integral', ['representation', 'integral_type', 'subdomain_id',
                                         'form_id', 'rank', 'geometric_dimension', 'topological_dimension',
                                         'entitytype', 'num_facets', 'num_vertices', 'needs_oriented',
                                         'enabled_coefficients', 'classnames', 'element_dimensions',
                                         'tensor_shape', 'quadrature_rules', 'coefficient_numbering',
                                         'coefficient_offsets', 'params', 'unique_tables', 'unique_table_types',
                                         'piecewise_ir', 'varying_irs', 'all_num_points', 'classname',
                                         'prefix', 'integrals_metadata', 'integral_metadata'])
ir_tabulate_dof_coordinates = namedtuple('ir_tabulate_dof_coordinates', ['tdim', 'gdim', 'points', 'cell_shape'])
ir_evaluate_dof = namedtuple('ir_evaluate_dof', ['mappings', 'reference_value_size', 'physical_value_size',
                                                 'geometric_dimension', 'topological_dimension', 'dofs',
                                                 'physical_offsets', 'cell_shape'])

ir_data = namedtuple('ir_data', ['elements', 'dofmaps', 'coordinate_mappings', 'integrals', 'forms'])


def make_finite_element_jit_classname(ufl_element, tag, parameters):
    assert isinstance(ufl_element, ufl.FiniteElementBase)
    sig = classname.compute_signature([ufl_element], tag, parameters)
    return classname.make_name("ffc_element_{}".format(sig), "finite_element", "main")


def make_dofmap_jit_classname(ufl_element, tag, parameters):
    assert isinstance(ufl_element, ufl.FiniteElementBase)
    sig = classname.compute_signature([ufl_element], tag, parameters)
    return classname.make_name("ffc_element_{}".format(sig), "dofmap", "main")


def make_coordinate_mapping_jit_classname(ufl_element, tag, parameters):
    assert isinstance(ufl_element, ufl.FiniteElementBase)
    sig = classname.compute_signature([ufl_element], tag, parameters, coordinate_mapping=True)
    return classname.make_name("ffc_coordinate_mapping_{}".format(sig), "coordinate_mapping", "main")


def make_all_element_classnames(prefix, elements, coordinate_elements, parameters):
    # Make unique classnames to match separately jit-compiled module
    classnames = {
        "finite_element": {e: make_finite_element_jit_classname(e, prefix, parameters)
                           for e in elements},
        "dofmap": {e: make_dofmap_jit_classname(e, prefix, parameters)
                   for e in elements},
        "coordinate_mapping":
        {e: make_coordinate_mapping_jit_classname(e, prefix, parameters)
         for e in coordinate_elements},
    }
    return classnames


def compute_ir(analysis: namedtuple, object_names, prefix, parameters):
    """Compute intermediate representation.

    """
    logger.info("Compiler stage 2: Computing intermediate representation")

    # Construct classnames for all element objects and coordinate mappings
    classnames = make_all_element_classnames(prefix, analysis.unique_elements,
                                             analysis.unique_coordinate_elements, parameters)

    # Compute representation of elements
    logger.info("Computing representation of {} elements".format(len(analysis.unique_elements)))
    ir_elements = [
        _compute_element_ir(e, analysis.element_numbers, classnames, parameters) for e in analysis.unique_elements
    ]

    # Compute representation of dofmaps
    logger.info("Computing representation of {} dofmaps".format(len(analysis.unique_elements)))
    ir_dofmaps = [
        _compute_dofmap_ir(e, analysis.element_numbers, classnames, parameters) for e in analysis.unique_elements
    ]

    # Compute representation of coordinate mappings
    logger.info("Computing representation of {} coordinate mappings".format(
        len(analysis.unique_coordinate_elements)))
    ir_coordinate_mappings = [
        _compute_coordinate_mapping_ir(e, analysis.element_numbers, classnames, parameters)
        for e in analysis.unique_coordinate_elements
    ]

    # Compute and flatten representation of integrals
    logger.info("Computing representation of integrals")
    irs = [
        _compute_integral_ir(fd, i, prefix, analysis.element_numbers, classnames, parameters)
        for (i, fd) in enumerate(analysis.form_data)
    ]
    ir_integrals = list(itertools.chain(*irs))

    # Compute representation of forms
    logger.info("Computing representation of forms")
    ir_forms = [
        _compute_form_ir(fd, i, prefix, analysis.element_numbers,
                         classnames, object_names, parameters)
        for (i, fd) in enumerate(analysis.form_data)
    ]

    return ir_data(elements=ir_elements, dofmaps=ir_dofmaps,
                   coordinate_mappings=ir_coordinate_mappings,
                   integrals=ir_integrals, forms=ir_forms)


def _compute_element_ir(ufl_element, element_numbers, classnames, parameters):
    """Compute intermediate representation of element."""
    # Create FIAT element
    fiat_element = create_element(ufl_element)
    cell = ufl_element.cell()
    cellname = cell.cellname()

    # Store id
    ir = {"id": element_numbers[ufl_element]}
    ir["classname"] = classnames["finite_element"][ufl_element]

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

    ir["evaluate_basis"] = _evaluate_basis(ufl_element, fiat_element, parameters["epsilon"])
    ir["evaluate_dof"] = _evaluate_dof(ufl_element, fiat_element)
    ir["tabulate_dof_coordinates"] = _tabulate_dof_coordinates(ufl_element, fiat_element)
    ir["num_sub_elements"] = ufl_element.num_sub_elements()
    ir["create_sub_element"] = [classnames["finite_element"][e] for e in ufl_element.sub_elements()]

    return ir_element(**ir)


def _compute_dofmap_ir(ufl_element, element_numbers, classnames, parameters):
    """Compute intermediate representation of dofmap."""
    # Create FIAT element
    fiat_element = create_element(ufl_element)

    # Precompute repeatedly used items
    num_dofs_per_entity = _num_dofs_per_entity(fiat_element)
    entity_dofs = fiat_element.entity_dofs()

    # Store id
    ir = {"id": element_numbers[ufl_element]}
    ir["classname"] = classnames["dofmap"][ufl_element]

    # Compute data for each function
    ir["signature"] = "FFC dofmap for " + repr(ufl_element)
    ir["num_global_support_dofs"] = _num_global_support_dofs(fiat_element)
    ir["num_element_support_dofs"] = fiat_element.space_dimension() - ir["num_global_support_dofs"]
    ir["num_entity_dofs"] = num_dofs_per_entity
    ir["tabulate_entity_dofs"] = (entity_dofs, num_dofs_per_entity)
    ir["num_sub_dofmaps"] = ufl_element.num_sub_elements()
    ir["create_sub_dofmap"] = [classnames["dofmap"][e] for e in ufl_element.sub_elements()]

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
                                   element_numbers,
                                   classnames,
                                   parameters):
    """Compute intermediate representation of coordinate mapping."""
    cell = ufl_coordinate_element.cell()
    cellname = cell.cellname()

    assert ufl_coordinate_element.value_shape() == (cell.geometric_dimension(), )

    # Compute element values via fiat element
    tables = _tabulate_coordinate_mapping_basis(ufl_coordinate_element)

    # Store id
    ir = {"id": element_numbers[ufl_coordinate_element]}
    ir["classname"] = classnames["coordinate_mapping"][ufl_coordinate_element]

    # Compute data for each function
    ir["signature"] = "FFC coordinate_mapping from " + repr(ufl_coordinate_element)
    ir["cell_shape"] = cellname
    ir["topological_dimension"] = cell.topological_dimension()
    ir["geometric_dimension"] = ufl_coordinate_element.value_size()

    ir["create_coordinate_finite_element"] = classnames["finite_element"][ufl_coordinate_element]
    ir["create_coordinate_dofmap"] = classnames["dofmap"][ufl_coordinate_element]

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

    # Get classnames for coordinate element and its scalar subelement:
    ir["coordinate_finite_element_classname"] = classnames["finite_element"][ufl_coordinate_element]
    ir["scalar_coordinate_finite_element_classname"] = classnames["finite_element"][
        ufl_coordinate_element.sub_elements()[0]]

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


def _compute_integral_ir(form_data, form_index, prefix, element_numbers, classnames, parameters):
    """Compute intermediate represention for form integrals."""
    if form_data.representation == "uflacs":
        from ffc.ir.uflacs.uflacsrepresentation import compute_integral_ir
    elif form_data.representation == "tsfc":
        from ffc.ir.tsfcrepresentation import compute_integral_ir
    else:
        raise RuntimeError("Unknown representation: {}".format(form_data.representation))

    # Iterate over integrals
    irs = []
    for itg_data in form_data.integral_data:
        # FIXME: Can we remove form_index?
        # Compute representation
        ir = compute_integral_ir(itg_data, form_data, form_index, element_numbers, classnames, parameters)

        # Build classname
        ir["classname"] = classname.make_integral_name(prefix, itg_data.integral_type, form_index,
                                                       itg_data.subdomain_id)
        ir["classnames"] = classnames  # FIXME XXX: Use this everywhere needed?

        # Storing prefix here for reconstruction of classnames on code
        # generation side
        ir["prefix"] = prefix  # FIXME: Drop this?

        # Store metadata for later reference (eg. printing as comment)
        # NOTE: We make a commitment not to modify it!
        ir["integrals_metadata"] = itg_data.metadata
        ir["integral_metadata"] = [integral.metadata() for integral in itg_data.integrals]

        irs.append(ir_integral(**ir))

    return irs


def _compute_form_ir(form_data, form_id, prefix, element_numbers,
                     classnames, object_names, parameters):
    """Compute intermediate representation of form."""

    # Store id
    ir = {"id": form_id}

    # Storing prefix here for reconstruction of classnames on code
    # generation side
    ir["prefix"] = prefix

    # Compute common data
    ir["classname"] = classname.make_name(prefix, "form", form_id)

    ir["signature"] = form_data.original_form.signature()

    ir["rank"] = len(form_data.original_form.arguments())
    ir["num_coefficients"] = len(form_data.reduced_coefficients)

    ir["coefficient_names"] = [object_names.get(id(obj), "w%d" % j)
                               for j, obj in enumerate(form_data.reduced_coefficients)]

    ir["original_coefficient_position"] = form_data.original_coefficient_positions

    # TODO: Remove create_coordinate_{finite_element,dofmap} and access
    # through coordinate_mapping instead in dolfin, when that's in place
    ir["create_coordinate_finite_element"] = [
        classnames["finite_element"][e] for e in form_data.coordinate_elements
    ]
    ir["create_coordinate_dofmap"] = [
        classnames["dofmap"][e] for e in form_data.coordinate_elements
    ]
    ir["create_coordinate_mapping"] = [
        classnames["coordinate_mapping"][e] for e in form_data.coordinate_elements
    ]
    ir["create_finite_element"] = [
        classnames["finite_element"][e]
        for e in form_data.argument_elements + form_data.coefficient_elements
    ]
    ir["create_dofmap"] = [
        classnames["dofmap"][e]
        for e in form_data.argument_elements + form_data.coefficient_elements
    ]

    # Create integral ids and names using form prefix (integrals are
    # always generated as part of form so don't get their own prefix)
    for integral_type in ufc_integral_types:
        irdata = _create_foo_integral(prefix, form_id, integral_type, form_data)
        ir["create_{}_integral".format(integral_type)] = irdata
        ir["get_{}_integral_ids".format(integral_type)] = irdata

    return ir_form(**ir)


def _generate_reference_offsets(fiat_element, offset=0):
    """Generate offsets: i.e value offset for each basis function
    relative to a reference element representation."""
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
    """Generate offsets: i.e value offset for each basis function
    relative to a physical element representation."""
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
    """Generate offsets: i.e value offset for each basis function
    relative to a physical element representation."""
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
        "needs_oriented": needs_oriented_jacobian(fiat_element),
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
        classnames += [classname.make_integral_name(prefix, integral_type, form_id, 'otherwise')]

    for itg_data in form_data.integral_data:
        if isinstance(itg_data.subdomain_id, int):
            if itg_data.subdomain_id < 0:
                raise ValueError("Integral subdomain ID must be non-negative, not {}".format(itg_data.subdomain_id))
            if (itg_data.integral_type == integral_type):
                subdomain_ids += [itg_data.subdomain_id]
                classnames += [classname.make_integral_name(prefix, integral_type,
                                                            form_id, itg_data.subdomain_id)]

    return subdomain_ids, classnames


def all_elements(fiat_element):
    if isinstance(fiat_element, MixedElement):
        return fiat_element.elements()
    return [fiat_element]


def _num_dofs_per_entity(fiat_element):
    """Compute list of integers representing the number of dofs
    associated with a single mesh entity.

    Example: Lagrange of degree 3 on triangle: [1, 2, 1]

    """
    entity_dofs = fiat_element.entity_dofs()
    return [len(entity_dofs[e][0]) for e in sorted(entity_dofs.keys())]


def uses_integral_moments(fiat_element):
    """True if element uses integral moments for its degrees of freedom.

    """
    integrals = set(["IntegralMoment", "FrobeniusIntegralMoment"])
    tags = set([L.get_type_tag() for L in fiat_element.dual_basis() if L])
    return len(integrals & tags) > 0


def needs_oriented_jacobian(fiat_element):
    # Check whether this element needs an oriented jacobian (only
    # contravariant piolas seem to need it)
    return "contravariant piola" in fiat_element.mapping()
