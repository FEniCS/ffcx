# -*- coding: utf-8 -*-
"""
Compiler stage 2: Code representation
-------------------------------------

This module computes intermediate representations of forms,
elements and dofmaps. For each UFC function, we extract the
data needed for code generation at a later stage.

The representation should conform strictly to the naming and
order of functions in UFC. Thus, for code generation of the
function "foo", one should only need to use the data stored
in the intermediate representation under the key "foo".
"""

# Copyright (C) 2009-2017 Anders Logg, Martin Sandve Alnæs, Marie E. Rognes,
# Kristian B. Oelgaard, and others
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.

# Python modules
from itertools import chain
import numpy

# Import UFL
import ufl

# FFC modules
from ffc.utils import compute_permutations, product
from ffc.log import info, error, begin, end
from ffc.fiatinterface import create_element, reference_cell
from ffc.fiatinterface import EnrichedElement, HDivTrace, MixedElement, SpaceOfReals, QuadratureElement
from ffc.classname import make_classname, make_integral_classname

# List of supported integral types
ufc_integral_types = ("cell",
                      "exterior_facet",
                      "interior_facet",
                      "vertex",
                      "custom",
                      "cutcell",
                      "interface",
                      "overlap")


def pick_representation(representation):
    "Return one of the specialized code generation modules from a representation string."
    if representation == "quadrature":
        from ffc import quadrature as r
    elif representation == "uflacs":
        from ffc import uflacs as r
    elif representation == "tsfc":
        from ffc import tsfc as r
    else:
        error("Unknown representation: %s" % str(representation))
    return r


def make_finite_element_jit_classname(ufl_element, parameters):
    from ffc.jitcompiler import compute_jit_prefix  # FIXME circular file dependency
    kind, prefix = compute_jit_prefix(ufl_element, parameters)
    return make_classname(prefix, "finite_element", "main")


def make_dofmap_jit_classname(ufl_element, parameters):
    from ffc.jitcompiler import compute_jit_prefix  # FIXME circular file dependency
    kind, prefix = compute_jit_prefix(ufl_element, parameters)
    return make_classname(prefix, "dofmap", "main")


def make_coordinate_mapping_jit_classname(ufl_mesh, parameters):
    from ffc.jitcompiler import compute_jit_prefix  # FIXME circular file dependency
    kind, prefix = compute_jit_prefix(ufl_mesh, parameters, kind="coordinate_mapping")
    return make_classname(prefix, "coordinate_mapping", "main")


def make_all_element_classnames(prefix, elements, coordinate_elements,
                                element_numbers, parameters, jit):
    if jit:
        # Make unique classnames to match separately jit-compiled
        # module
        classnames = {
            "finite_element": {
                e: make_finite_element_jit_classname(e, parameters)
                for e in elements },
            "dofmap": {
                e: make_dofmap_jit_classname(e, parameters)
                for e in elements },
            "coordinate_mapping": {
                e: make_coordinate_mapping_jit_classname(e, parameters)
                for e in coordinate_elements },
            }
    else:
        # Make unique classnames only within this module (using a
        # shared prefix and element numbers that are only unique
        # within this module)
        classnames = {
            "finite_element": {
                e: make_classname(prefix, "finite_element", element_numbers[e])
                for e in elements },
            "dofmap": {
                e: make_classname(prefix, "dofmap", element_numbers[e])
                for e in elements },
            "coordinate_mapping": {
                e: make_classname(prefix, "coordinate_mapping", element_numbers[e])
                for e in coordinate_elements },
            }
    return classnames


def compute_ir(analysis, prefix, parameters, jit=False):
    "Compute intermediate representation."

    begin("Compiler stage 2: Computing intermediate representation")

    # Set code generation parameters (this is not actually a 'formatting'
    # parameter, used for table value clamping as well)
    # FIXME: Global state?!
    #    set_float_formatting(parameters["precision"])

    # Extract data from analysis
    form_datas, elements, element_numbers, coordinate_elements = analysis

    # Construct classnames for all element objects and coordinate mappings
    classnames = make_all_element_classnames(prefix, elements,
                                             coordinate_elements,
                                             element_numbers,
                                             parameters, jit)

    # Skip processing elements if jitting forms
    # NB! it's important that this happens _after_ the element numbers and classnames
    # above have been created.
    if jit and form_datas:
        # While we may get multiple forms during command line action,
        # not so during jit
        assert len(form_datas) == 1, "Expecting only one form data instance during jit."
        # Drop some processing
        elements = []
        coordinate_elements = []
    elif jit and coordinate_elements:
        # While we may get multiple coordinate elements during command
        # line action, or during form jit, not so during coordinate
        # mapping jit
        assert len(coordinate_elements) == 1, "Expecting only one form data instance during jit."
        # Drop some processing
        elements = []
    elif jit and elements:
        # Assuming a topological sorting of the elements,
        # only process the last (main) element from here on
        elements = [elements[-1]]

    # Compute representation of elements
    info("Computing representation of %d elements" % len(elements))
    ir_elements = [_compute_element_ir(e, element_numbers, classnames, parameters, jit)
                   for e in elements]

    # Compute representation of dofmaps
    info("Computing representation of %d dofmaps" % len(elements))
    ir_dofmaps = [_compute_dofmap_ir(e, element_numbers, classnames, parameters, jit)
                  for e in elements]

    # Compute representation of coordinate mappings
    info("Computing representation of %d coordinate mappings" % len(coordinate_elements))
    ir_coordinate_mappings = [_compute_coordinate_mapping_ir(e, element_numbers, classnames, parameters, jit)
                              for e in coordinate_elements]

    # Compute and flatten representation of integrals
    info("Computing representation of integrals")
    irs = [_compute_integral_ir(fd, form_id, prefix, element_numbers, classnames, parameters, jit)
           for (form_id, fd) in enumerate(form_datas)]
    ir_integrals = list(chain(*irs))

    # Compute representation of forms
    info("Computing representation of forms")
    ir_forms = [_compute_form_ir(fd, form_id, prefix, element_numbers, classnames, parameters, jit)
                for (form_id, fd) in enumerate(form_datas)]

    end()

    return ir_elements, ir_dofmaps, ir_coordinate_mappings, ir_integrals, ir_forms


def _compute_element_ir(ufl_element, element_numbers, classnames, parameters, jit):
    "Compute intermediate representation of element."

    # Create FIAT element
    fiat_element = create_element(ufl_element)
    cell = ufl_element.cell()
    cellname = cell.cellname()

    # Store id
    ir = {"id": element_numbers[ufl_element]}
    ir["classname"] = classnames["finite_element"][ufl_element]

    # Remember jit status
    ir["jit"] = jit

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
    ir["interpolate_vertex_values"] = _interpolate_vertex_values(ufl_element,
                                                                 fiat_element)
    ir["tabulate_dof_coordinates"] = _tabulate_dof_coordinates(ufl_element,
                                                               fiat_element)
    ir["num_sub_elements"] = ufl_element.num_sub_elements()
    ir["create_sub_element"] = [classnames["finite_element"][e]
                                for e in ufl_element.sub_elements()]

    # debug_ir(ir, "finite_element")

    return ir


def _compute_dofmap_ir(ufl_element, element_numbers, classnames, parameters, jit=False):
    "Compute intermediate representation of dofmap."

    # Create FIAT element
    fiat_element = create_element(ufl_element)
    cell = ufl_element.cell()

    # Precompute repeatedly used items
    num_dofs_per_entity = _num_dofs_per_entity(fiat_element)
    entity_dofs = fiat_element.entity_dofs()
    facet_dofs = _tabulate_facet_dofs(fiat_element, cell)
    entity_closure_dofs, num_dofs_per_entity_closure = \
        _tabulate_entity_closure_dofs(fiat_element, cell)

    # Store id
    ir = {"id": element_numbers[ufl_element]}
    ir["classname"] = classnames["dofmap"][ufl_element]

    # Remember jit status
    ir["jit"] = jit

    # Compute data for each function
    ir["signature"] = "FFC dofmap for " + repr(ufl_element)
    ir["needs_mesh_entities"] = _needs_mesh_entities(fiat_element)
    ir["topological_dimension"] = cell.topological_dimension()
    ir["geometric_dimension"] = cell.geometric_dimension()
    ir["global_dimension"] = _global_dimension(fiat_element)
    ir["num_global_support_dofs"] = _num_global_support_dofs(fiat_element)
    ir["num_element_support_dofs"] = fiat_element.space_dimension() - ir["num_global_support_dofs"]
    ir["num_element_dofs"] = fiat_element.space_dimension()
    ir["num_facet_dofs"] = len(facet_dofs[0])
    ir["num_entity_dofs"] = num_dofs_per_entity
    ir["num_entity_closure_dofs"] = num_dofs_per_entity_closure
    ir["tabulate_dofs"] = _tabulate_dofs(fiat_element, cell)
    ir["tabulate_facet_dofs"] = facet_dofs
    ir["tabulate_entity_dofs"] = (entity_dofs, num_dofs_per_entity)
    ir["tabulate_entity_closure_dofs"] = (entity_closure_dofs, entity_dofs, num_dofs_per_entity)
    ir["num_sub_dofmaps"] = ufl_element.num_sub_elements()
    ir["create_sub_dofmap"] = [classnames["dofmap"][e]
                               for e in ufl_element.sub_elements()]

    return ir


_midpoints = {
    "interval": (0.5,),
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
    origo = (0.0,) * tdim
    midpoint = cell_midpoint(cell)

    # Tabulate basis
    t0 = fiat_element.tabulate(1, [origo])
    tm = fiat_element.tabulate(1, [midpoint])

    # Get basis values at cell origo
    tables["x0"] = t0[(0,) * tdim][:, 0]

    # Get basis values at cell midpoint
    tables["xm"] = tm[(0,) * tdim][:, 0]

    # Single direction derivatives, e.g. [(1,0), (0,1)] in 2d
    derivatives = [(0,) * i + (1,) + (0,) * (tdim - 1 - i) for i in range(tdim)]

    # Get basis derivative values at cell origo
    tables["J0"] = numpy.asarray([t0[d][:, 0] for d in derivatives])

    # Get basis derivative values at cell midpoint
    tables["Jm"] = numpy.asarray([tm[d][:, 0] for d in derivatives])

    return tables


def _compute_coordinate_mapping_ir(ufl_coordinate_element, element_numbers,
                                   classnames, parameters, jit=False):
    "Compute intermediate representation of coordinate mapping."

    cell = ufl_coordinate_element.cell()
    cellname = cell.cellname()

    assert ufl_coordinate_element.value_shape() == (cell.geometric_dimension(),)

    # Compute element values via fiat element
    tables = _tabulate_coordinate_mapping_basis(ufl_coordinate_element)

    # Store id
    ir = {"id": element_numbers[ufl_coordinate_element]}
    ir["classname"] = classnames["coordinate_mapping"][ufl_coordinate_element]

    # Remember jit status
    ir["jit"] = jit

    # Compute data for each function
    ir["signature"] = "FFC coordinate_mapping from " + repr(ufl_coordinate_element)
    ir["cell_shape"] = cellname
    ir["topological_dimension"] = cell.topological_dimension()
    ir["geometric_dimension"] = ufl_coordinate_element.value_size()

    ir["create_coordinate_finite_element"] = \
        classnames["finite_element"][ufl_coordinate_element]
    ir["create_coordinate_dofmap"] = \
        classnames["dofmap"][ufl_coordinate_element]

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
    ir["coordinate_finite_element_classname"] = \
        classnames["finite_element"][ufl_coordinate_element]
    ir["scalar_coordinate_finite_element_classname"] = \
        classnames["finite_element"][ufl_coordinate_element.sub_elements()[0]]

    return ir


def _num_global_support_dofs(fiat_element):
    "Compute number of global support dofs."
    if not isinstance(fiat_element, MixedElement):
        if isinstance(fiat_element, SpaceOfReals):
            return 1
        return 0
    num_reals = 0
    for e in fiat_element.elements():
        if isinstance(e, SpaceOfReals):
            num_reals += 1
    return num_reals


def _global_dimension(fiat_element):
    "Compute intermediate representation for global_dimension."

    if not isinstance(fiat_element, MixedElement):
        if isinstance(fiat_element, SpaceOfReals):
            return ([], 1)
        return (_num_dofs_per_entity(fiat_element), 0)

    elements = []
    num_reals = 0
    for e in fiat_element.elements():
        if not isinstance(e, SpaceOfReals):
            elements += [e]
        else:
            num_reals += 1
    fiat_cell, = set(e.get_reference_element()
                     for e in fiat_element.elements())
    fiat_element = MixedElement(elements, ref_el=fiat_cell)
    return (_num_dofs_per_entity(fiat_element), num_reals)


def _needs_mesh_entities(fiat_element):
    "Compute intermediate representation for needs_mesh_entities."

    # Note: The dof map for Real elements does not depend on the mesh

    num_dofs_per_entity = _num_dofs_per_entity(fiat_element)
    if isinstance(fiat_element, SpaceOfReals):
        return [False for d in num_dofs_per_entity]
    else:
        return [d > 0 for d in num_dofs_per_entity]


def _compute_integral_ir(form_data, form_id, prefix, element_numbers, classnames,
                         parameters, jit):
    "Compute intermediate represention for form integrals."

    # For consistency, all jit objects now have classnames with postfix "main"
    if jit:
        assert form_id == 0
        form_id = "main"

    irs = []

    # Iterate over integrals
    for itg_data in form_data.integral_data:

        # Select representation
        # TODO: Is it possible to detach this metadata from
        # IntegralData? It's a bit strange from the ufl side.
        r = pick_representation(itg_data.metadata["representation"])

        # Compute representation
        ir = r.compute_integral_ir(itg_data,
                                   form_data,
                                   form_id,     # FIXME: Can we remove this?
                                   element_numbers,
                                   classnames,
                                   parameters)

        # Remember jit status
        ir["jit"] = jit

        # Build classname
        ir["classname"] = make_integral_classname(prefix, itg_data.integral_type,
                                                  form_id, itg_data.subdomain_id)

        ir["classnames"] = classnames  # FIXME XXX: Use this everywhere needed?

        # Storing prefix here for reconstruction of classnames on code
        # generation side
        ir["prefix"] = prefix  # FIXME: Drop this?

        # Store metadata for later reference (eg. printing as comment)
        # NOTE: We make a commitment not to modify it!
        ir["integrals_metadata"] = itg_data.metadata
        ir["integral_metadata"] = [integral.metadata()
                                   for integral in itg_data.integrals]

        # Append representation
        irs.append(ir)

    return irs


def _compute_form_ir(form_data, form_id, prefix, element_numbers,
                     classnames, parameters, jit=False):
    "Compute intermediate representation of form."

    # For consistency, all jit objects now have classnames with postfix "main"
    if jit:
        assert form_id == 0
        form_id = "main"

    # Store id
    ir = {"id": form_id}

    # Storing prefix here for reconstruction of classnames on code
    # generation side
    ir["prefix"] = prefix

    # Remember jit status
    ir["jit"] = jit

    # Compute common data
    ir["classname"] = make_classname(prefix, "form", form_id)

    # ir["members"] = None # unused
    # ir["constructor"] = None # unused
    # ir["destructor"] = None # unused
    ir["signature"] = form_data.original_form.signature()

    ir["rank"] = len(form_data.original_form.arguments())
    ir["num_coefficients"] = len(form_data.reduced_coefficients)
    ir["original_coefficient_position"] = form_data.original_coefficient_positions

    # TODO: Remove create_coordinate_{finite_element,dofmap} and
    # access through coordinate_mapping instead in dolfin, when that's
    # in place
    ir["create_coordinate_finite_element"] = [
        classnames["finite_element"][e]
        for e in form_data.coordinate_elements
        ]
    ir["create_coordinate_dofmap"] = [
        classnames["dofmap"][e]
        for e in form_data.coordinate_elements
        ]
    ir["create_coordinate_mapping"] = [
        classnames["coordinate_mapping"][e]
        for e in form_data.coordinate_elements
        ]
    ir["create_finite_element"] = [
        classnames["finite_element"][e]
        for e in form_data.argument_elements + form_data.coefficient_elements
        ]
    ir["create_dofmap"] = [
        classnames["dofmap"][e]
        for e in form_data.argument_elements + form_data.coefficient_elements
        ]

    # Create integral ids and names using form prefix
    # (integrals are always generated as part of form so don't get
    # their own prefix)
    for integral_type in ufc_integral_types:
        ir["max_%s_subdomain_id" % integral_type] = \
            form_data.max_subdomain_ids.get(integral_type, 0)
        ir["has_%s_integrals" % integral_type] = \
            _has_foo_integrals(prefix, form_id, integral_type, form_data)
        ir["create_%s_integral" % integral_type] = \
            _create_foo_integral(prefix, form_id, integral_type, form_data)
        ir["create_default_%s_integral" % integral_type] = \
            _create_default_foo_integral(prefix, form_id, integral_type, form_data)

    return ir


# --- Computation of intermediate representation for non-trivial functions ---

def _generate_reference_offsets(fiat_element, offset=0):
    """Generate offsets: i.e value offset for each basis function
    relative to a reference element representation."""

    if isinstance(fiat_element, MixedElement):
        offsets = []
        for e in fiat_element.elements():
            offsets += _generate_reference_offsets(e, offset)
            # NB! This is the fiat element and therefore value_shape
            # means reference_value_shape
            offset += product(e.value_shape())
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

    # Refer to reference if gdim == tdim. This is a hack to support
    # more stuff (in particular restricted elements)
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
    "Compute intermediate representation of evaluate_dof."
    cell = ufl_element.cell()
    if fiat_element.is_nodal():
        dofs = [L.pt_dict for L in fiat_element.dual_basis()]
    else:
        dofs = [None] * fiat_element.space_dimension()
    return {"mappings": fiat_element.mapping(),
            "reference_value_size": ufl_element.reference_value_size(),
            "physical_value_size": ufl_element.value_size(),
            "geometric_dimension": cell.geometric_dimension(),
            "topological_dimension": cell.topological_dimension(),
            "dofs": dofs,
            "physical_offsets": _generate_physical_offsets(ufl_element),
            "cell_shape": cell.cellname()}


def _extract_elements(fiat_element):
    new_elements = []
    if isinstance(fiat_element, (MixedElement, EnrichedElement)):
        for e in fiat_element.elements():
            new_elements += _extract_elements(e)
    else:
        new_elements.append(fiat_element)
    return new_elements


def _evaluate_basis(ufl_element, fiat_element, epsilon):
    "Compute intermediate representation for evaluate_basis."
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

    # Handle QuadratureElement, not supported because the basis is
    # only defined at the dof coordinates where the value is 1, so not
    # very interesting.
    for e in elements:
        if isinstance(e, QuadratureElement):
            return "Function not supported/implemented for QuadratureElement."
        if isinstance(e, HDivTrace):
            return "Function not supported for Trace elements"

    # Skip this function for TensorProductElement if get_coeffs is not implemented
    for e in elements:
        try:
            e.get_coeffs()
        except NotImplementedError:
            return "Function is not supported/implemented."

    # Initialise data with 'global' values.
    data = {"reference_value_size": ufl_element.reference_value_size(),
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
        num_components = product(e.value_shape())
        coeffs = e.get_coeffs()
        num_expansion_members = e.get_num_members(e.degree())
        dmats = e.dmats()

        # Clamp dmats zeros
        dmats = numpy.asarray(dmats)
        dmats[numpy.where(numpy.isclose(dmats, 0.0, rtol=epsilon, atol=epsilon))] = 0.0

        # Extracted parts of dd below that are common for the element
        # here.  These dict entries are added to each dof_data dict
        # for each dof, because that's what the code generation
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
                coefficients = [coeffs[i][c]
                                for c in range(num_components)]
            elif value_rank == 2:
                # Handle coefficients for tensor valued basis elements.
                # [Regge]
                coefficients = [coeffs[i][p][q]
                                for p in range(e.value_shape()[0])
                                for q in range(e.value_shape()[1])]
            else:
                error("Unknown situation with num_components > 1")

            # Clamp coefficient zeros
            coefficients = numpy.asarray(coefficients)
            coefficients[numpy.where(numpy.isclose(coefficients, 0.0, rtol=epsilon, atol=epsilon))] = 0.0

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
    "Compute intermediate representation of tabulate_dof_coordinates."
    if uses_integral_moments(element):
        return {}

    # Bail out if any dual basis member is missing (element is not nodal),
    # this is strictly not necessary but simpler
    if any(L is None for L in element.dual_basis()):
        return {}

    cell = ufl_element.cell()

    data = {}
    data["tdim"] = cell.topological_dimension()
    data["gdim"] = cell.geometric_dimension()
    data["points"] = [sorted(L.pt_dict.keys())[0] for L in element.dual_basis()]
    data["cell_shape"] = cell.cellname()
    return data


def _tabulate_dofs(element, cell):
    "Compute intermediate representation of tabulate_dofs."

    if isinstance(element, SpaceOfReals):
        return None

    # Extract number of dofs per entity for each element
    elements = all_elements(element)
    num_dofs_per_element = [_num_dofs_per_entity(e) for e in elements]

    # Extract local dof numbers per entity for each element
    all_entity_dofs = [e.entity_dofs() for e in elements]
    dofs_per_element = [[[list(dofs[dim][entity])
                          for entity in sorted(dofs[dim].keys())]
                         for dim in sorted(dofs.keys())]
                        for dofs in all_entity_dofs]

    # Check whether we need offset
    multiple_entities = any([sum(n > 0 for n in num_dofs) - 1
                             for num_dofs in num_dofs_per_element])
    need_offset = len(elements) > 1 or multiple_entities

    num_dofs_per_element = [e.space_dimension() for e in elements]

    # Handle global "elements"
    fakes = [isinstance(e, SpaceOfReals) for e in elements]

    return (dofs_per_element, num_dofs_per_element, need_offset, fakes)


def _tabulate_facet_dofs(element, cell):
    "Compute intermediate representation of tabulate_facet_dofs."

    # Get topological dimension
    D = cell.topological_dimension()

    # Get the number of facets
    num_facets = cell.num_facets()

    # Make list of dofs
    facet_dofs = list(element.entity_closure_dofs()[D-1].values())

    assert num_facets == len(facet_dofs)

    facet_dofs = [facet_dofs[facet] for facet in range(num_facets)]

    return facet_dofs


def _tabulate_entity_closure_dofs(element, cell):
    "Compute intermediate representation of tabulate_entity_closure_dofs."

    # Get topological dimension
    D = cell.topological_dimension()

    # Get entity closure dofs from FIAT element
    fiat_entity_closure_dofs = element.entity_closure_dofs()

    entity_closure_dofs = {}
    for d0 in sorted(fiat_entity_closure_dofs.keys()):
        for e0 in sorted(fiat_entity_closure_dofs[d0].keys()):
            entity_closure_dofs[(d0, e0)] = fiat_entity_closure_dofs[d0][e0]

    num_entity_closure_dofs = [len(fiat_entity_closure_dofs[d0][0]) for d0 in sorted(fiat_entity_closure_dofs.keys())]

    return entity_closure_dofs, num_entity_closure_dofs


def _interpolate_vertex_values(ufl_element, fiat_element):
    "Compute intermediate representation of interpolate_vertex_values."

    # Check for QuadratureElement
    for e in all_elements(fiat_element):
        if isinstance(e, QuadratureElement):
            return "Function is not supported/implemented for QuadratureElement."
        if isinstance(e, HDivTrace):
            return "Function is not implemented for HDivTrace."

    cell = ufl_element.cell()
    cellname = cell.cellname()
    tdim = cell.topological_dimension()
    gdim = cell.geometric_dimension()

    ir = {}
    ir["geometric_dimension"] = gdim
    ir["topological_dimension"] = tdim

    # Check whether computing the Jacobian is necessary
    mappings = fiat_element.mapping()
    ir["needs_jacobian"] = any("piola" in m for m in mappings)
    ir["needs_oriented"] = needs_oriented_jacobian(fiat_element)

    # See note in _evaluate_dofs
    ir["reference_value_size"] = ufl_element.reference_value_size()
    ir["physical_value_size"] = ufl_element.value_size()

    # Get vertices of reference cell
    fiat_cell = reference_cell(cellname)
    vertices = fiat_cell.get_vertices()

    # Compute data for each constituent element
    all_fiat_elm = all_elements(fiat_element)
    ir["element_data"] = [
        {
            # NB! value_shape of fiat element e means reference_value_shape
           "reference_value_size": product(e.value_shape()),

           # FIXME: THIS IS A BUG:
           "physical_value_size": product(e.value_shape()),  # FIXME: Get from corresponding ufl element?

           "basis_values": e.tabulate(0, vertices)[(0,) * tdim].transpose(),
           "mapping": e.mapping()[0],
           "space_dim": e.space_dimension(),
        }
        for e in all_fiat_elm]

    # FIXME: Temporary hack!
    if len(ir["element_data"]) == 1:
        ir["element_data"][0]["physical_value_size"] = ir["physical_value_size"]

    # Consistency check, related to note in _evaluate_dofs
    # This will fail for e.g. (RT1 x DG0) on a manifold because of the above bug
    if sum(data["physical_value_size"] for data in ir["element_data"]) != ir["physical_value_size"]:
        ir = "Failed to set physical value size correctly for subelements."
    elif sum(data["reference_value_size"] for data in ir["element_data"]) != ir["reference_value_size"]:
        ir = "Failed to set reference value size correctly for subelements."

    return ir


def _has_foo_integrals(prefix, form_id, integral_type, form_data):
    "Compute intermediate representation of has_foo_integrals."
    v = (form_data.max_subdomain_ids.get(integral_type, 0) > 0
         or _create_default_foo_integral(prefix, form_id, integral_type, form_data) is not None)
    return bool(v)


def _create_foo_integral(prefix, form_id, integral_type, form_data):
    "Compute intermediate representation of create_foo_integral."
    subdomain_ids = [itg_data.subdomain_id
                     for itg_data in form_data.integral_data
                     if (itg_data.integral_type == integral_type
                         and isinstance(itg_data.subdomain_id, int))]
    classnames = [make_integral_classname(prefix, integral_type, form_id, subdomain_id)
                  for subdomain_id in subdomain_ids]
    return subdomain_ids, classnames


def _create_default_foo_integral(prefix, form_id, integral_type, form_data):
    "Compute intermediate representation of create_default_foo_integral."
    itg_data = [itg_data for itg_data in form_data.integral_data
                if (itg_data.integral_type == integral_type and
                    itg_data.subdomain_id == "otherwise")]
    if len(itg_data) > 1:
        error("Expecting at most one default integral of each type.")
    if itg_data:
        classname = make_integral_classname(prefix, integral_type, form_id, "otherwise")
        return classname
    else:
        return None


#--- Utility functions ---


def all_elements(fiat_element):
    if isinstance(fiat_element, MixedElement):
        return fiat_element.elements()
    return [fiat_element]


def _num_dofs_per_entity(fiat_element):
    """
    Compute list of integers representing the number of dofs
    associated with a single mesh entity.

    Example: Lagrange of degree 3 on triangle: [1, 2, 1]
    """
    entity_dofs = fiat_element.entity_dofs()
    return [len(entity_dofs[e][0]) for e in sorted(entity_dofs.keys())]


def uses_integral_moments(fiat_element):
    "True if element uses integral moments for its degrees of freedom."
    integrals = set(["IntegralMoment", "FrobeniusIntegralMoment"])
    tags = set([L.get_type_tag() for L in fiat_element.dual_basis() if L])
    return len(integrals & tags) > 0


def needs_oriented_jacobian(fiat_element):
    # Check whether this element needs an oriented jacobian
    # (only contravariant piolas seem to need it)
    return "contravariant piola" in fiat_element.mapping()
