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

# Copyright (C) 2009-2013 Anders Logg
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
#
# Modified by Marie E. Rognes 2010
# Modified by Kristian B. Oelgaard 2010
# Modified by Martin Alnaes, 2013
#
# First added:  2009-12-16
# Last changed: 2013-02-10

# Python modules
from itertools import chain

# Import UFL
import ufl
from ufl.classes import Measure

# FFC modules
from ffc.utils import compute_permutations, product
from ffc.log import info, error, begin, end, debug_ir, ffc_assert, warning
from ffc.fiatinterface import create_element, cellname_to_num_entities, reference_cell
from ffc.mixedelement import MixedElement
from ffc.enrichedelement import EnrichedElement, SpaceOfReals
from ffc.quadratureelement import QuadratureElement
from ffc.cpp import set_float_formatting

def pick_representation(representation):
    "Return one of the specialized code generation modules from a representation string."
    if representation == "quadrature":
        from ffc import quadrature
        r = quadrature
    elif representation == "tensor":
        from ffc import tensor
        r = tensor
    elif representation == "uflacs":
        from ffc import uflacsrepr
        r = uflacsrepr
    else:
        error("Unknown representation: %s" % str(representation))
    return r

not_implemented = None

def compute_ir(analysis, parameters):
    "Compute intermediate representation."

    begin("Compiler stage 2: Computing intermediate representation")

    # Set code generation parameters
    set_float_formatting(int(parameters["precision"]))

    # Extract data from analysis
    form_datas, elements, element_numbers = analysis

    # Compute representation of elements
    info("Computing representation of %d elements" % len(elements))
    ir_elements = [_compute_element_ir(e, i, element_numbers) \
                       for (i, e) in enumerate(elements)]

    # Compute representation of dofmaps
    info("Computing representation of %d dofmaps" % len(elements))
    ir_dofmaps = [_compute_dofmap_ir(e, i, element_numbers)
                      for (i, e) in enumerate(elements)]

    # Compute and flatten representation of integrals
    info("Computing representation of integrals")
    irs = [_compute_integral_ir(fd, i, parameters) \
               for (i, fd) in enumerate(form_datas)]
    ir_integrals = [ir for ir in chain(*irs) if not ir is None]

    # Compute representation of forms
    info("Computing representation of forms")
    ir_forms = [_compute_form_ir(fd, i, element_numbers) \
                    for (i, fd) in enumerate(form_datas)]

    end()

    return ir_elements, ir_dofmaps, ir_integrals, ir_forms

def _compute_element_ir(ufl_element, element_id, element_numbers):
    "Compute intermediate representation of element."

    # Create FIAT element
    element = create_element(ufl_element)
    cell = ufl_element.cell()

    # Store id
    ir = {"id": element_id}

    # Compute data for each function
    ir["signature"] = repr(ufl_element)
    ir["cell_shape"] = cell.cellname()
    ir["topological_dimension"] = cell.topological_dimension()
    ir["geometric_dimension"] = cell.geometric_dimension()
    ir["space_dimension"] = element.space_dimension()
    ir["value_rank"] = len(ufl_element.value_shape())
    ir["value_dimension"] = ufl_element.value_shape()
    ir["evaluate_basis"] = _evaluate_basis(ufl_element, element, cell)
    ir["evaluate_dof"] = _evaluate_dof(ufl_element, element, cell)
    ir["interpolate_vertex_values"] = _interpolate_vertex_values(ufl_element, element, cell)
    ir["num_sub_elements"] = ufl_element.num_sub_elements()
    ir["create_sub_element"] = _create_sub_foo(ufl_element, element_numbers)

    #debug_ir(ir, "finite_element")

    return ir

def _compute_dofmap_ir(ufl_element, element_id, element_numbers):
    "Compute intermediate representation of dofmap."

    # Create FIAT element
    element = create_element(ufl_element)
    cell = ufl_element.cell()

    # Precompute repeatedly used items
    num_dofs_per_entity = _num_dofs_per_entity(element)
    facet_dofs = _tabulate_facet_dofs(element, cell)

    # Store id
    ir = {"id": element_id}

    # Compute data for each function
    ir["signature"] = "FFC dofmap for " + repr(ufl_element)
    ir["needs_mesh_entities"] = _needs_mesh_entities(element)
    ir["topological_dimension"] = cell.topological_dimension()
    ir["geometric_dimension"] = cell.geometric_dimension()
    ir["global_dimension"] = _global_dimension(element)
    ir["local_dimension"] = element.space_dimension()
    ir["num_facet_dofs"] = len(facet_dofs[0])
    ir["num_entity_dofs"] = num_dofs_per_entity
    ir["tabulate_dofs"] = _tabulate_dofs(element, cell)
    ir["tabulate_facet_dofs"] = facet_dofs
    ir["tabulate_entity_dofs"] = (element.entity_dofs(), num_dofs_per_entity)
    ir["tabulate_coordinates"] = _tabulate_coordinates(ufl_element, element)
    ir["num_sub_dofmaps"] = ufl_element.num_sub_elements()
    ir["create_sub_dofmap"] = _create_sub_foo(ufl_element, element_numbers)

    #debug_ir(ir, "dofmap")

    return ir

def _global_dimension(element):
    "Compute intermediate representation for global_dimension."

    if not isinstance(element, MixedElement):
        if isinstance(element, SpaceOfReals):
            return ([], 1)
        return (_num_dofs_per_entity(element), 0)

    elements = []
    reals = []
    num_reals = 0
    for (i, e) in enumerate(element.elements()):
        if not isinstance(e, SpaceOfReals):
            elements += [e]
        else:
            num_reals += 1
    element = MixedElement(elements)
    return (_num_dofs_per_entity(element), num_reals)

def _needs_mesh_entities(element):
    "Compute intermediate representation for needs_mesh_entities."

    # Note: The dof map for Real elements does not depend on the mesh

    num_dofs_per_entity = _num_dofs_per_entity(element)
    if isinstance(element, SpaceOfReals):
        return [False for d in num_dofs_per_entity]
    else:
        return [d > 0 for d in num_dofs_per_entity]

def _compute_integral_ir(form_data, form_id, parameters):
    "Compute intermediate represention for form integrals."

    irs = []

    # Iterate over integrals
    for itg_data in form_data.integral_data:

        # Select representation
        # TODO: Is it possible to detach this metadata from IntegralData? It's a bit strange from the ufl side.
        r = pick_representation(itg_data.metadata["representation"])

        # Compute representation
        ir = r.compute_integral_ir(itg_data,
                                   form_data,
                                   form_id,
                                   parameters)

        # Append representation
        irs.append(ir)

    return irs

def _compute_form_ir(form_data, form_id, element_numbers):
    "Compute intermediate representation of form."

    # Store id
    ir = {"id": form_id}

    # Compute common data
    ir["classname"] = "FooForm"
    ir["members"] = not_implemented
    ir["constructor"] = not_implemented
    ir["destructor"] = not_implemented
    ir["signature"] = form_data.signature
    ir["rank"] = form_data.rank
    ir["num_coefficients"] = form_data.num_coefficients
    ir["num_cell_domains"] = form_data.num_sub_domains.get("cell",0)
    ir["num_exterior_facet_domains"] = form_data.num_sub_domains.get("exterior_facet",0)
    ir["num_interior_facet_domains"] = form_data.num_sub_domains.get("interior_facet",0)
    ir["num_point_domains"] = form_data.num_sub_domains.get("point",0)
    ir["has_cell_integrals"] = _has_foo_integrals("cell", form_data)
    ir["has_exterior_facet_integrals"] = _has_foo_integrals("exterior_facet", form_data)
    ir["has_interior_facet_integrals"] = _has_foo_integrals("interior_facet", form_data)
    ir["has_point_integrals"] = _has_foo_integrals("point", form_data)
    ir["create_finite_element"] = [element_numbers[e] for e in form_data.elements]
    ir["create_dofmap"] = [element_numbers[e] for e in form_data.elements]
    ir["create_cell_integral"] = _create_foo_integral("cell", form_data)
    ir["create_exterior_facet_integral"] = _create_foo_integral("exterior_facet", form_data)
    ir["create_interior_facet_integral"] = _create_foo_integral("interior_facet", form_data)
    ir["create_point_integral"] = _create_foo_integral("point", form_data)
    ir["create_default_cell_integral"] = _create_default_foo_integral("cell", form_data)
    ir["create_default_exterior_facet_integral"] = _create_default_foo_integral("exterior_facet", form_data)
    ir["create_default_interior_facet_integral"] = _create_default_foo_integral("interior_facet", form_data)
    ir["create_default_point_integral"] = _create_default_foo_integral("point", form_data)

    return ir

#--- Computation of intermediate representation for non-trivial functions ---

# FIXME: Move to FiniteElement/MixedElement
def _value_size(element):
    """Compute value size of element, aka the number of components.

    The value size of a scalar field is 1, the value size of a vector
    field (is the number of components), the value size of a higher
    dimensional tensor field is the product of the value_shape of the
    field. Recall that all mixed elements are flattened.
    """
    shape = element.value_shape()
    if shape == ():
        return 1
    return product(shape)

def _generate_reference_offsets(element, offset=0):
    """Generate offsets: i.e value offset for each basis function
    relative to a reference element representation."""
    offsets = []
    if isinstance(element, MixedElement):
        for e in element.elements():
            offsets += _generate_reference_offsets(e, offset)
            offset += _value_size(e)
    elif isinstance(element, EnrichedElement):
        for e in element.elements():
            offsets += _generate_reference_offsets(e, offset)
    else:
        offsets = [offset]*element.space_dimension()
    return offsets

def _generate_physical_offsets(ufl_element, offset=0):
    """Generate offsets: i.e value offset for each basis function
    relative to a physical element representation."""
    offsets = []

    # Refer to reference if gdim == tdim. This is a hack to support
    # more stuff (in particular restricted elements)
    gdim = ufl_element.cell().geometric_dimension()
    tdim = ufl_element.cell().topological_dimension()
    if (gdim == tdim):
        return _generate_reference_offsets(create_element(ufl_element))

    if isinstance(ufl_element, ufl.MixedElement):
        for e in ufl_element.sub_elements():
            offsets += _generate_physical_offsets(e, offset)
            offset += _value_size(e)
    elif isinstance(ufl_element, ufl.EnrichedElement):
        for e in ufl_element._elements:
            offsets += _generate_physical_offsets(e, offset)
    elif isinstance(ufl_element, ufl.FiniteElement):
        element = create_element(ufl_element)
        offsets = [offset]*element.space_dimension()
    else:
        raise NotImplementedError, \
            "This element combination is not implemented"
    return offsets

def _evaluate_dof(ufl_element, element, cell):
    "Compute intermediate representation of evaluate_dof."

    # With regard to reference_value_size vs physical_value_size: Note
    # that 'element' is the FFC/FIAT representation of the finite
    # element, while 'ufl_element' is the UFL representation. In
    # particular, UFL only knows about physical dimensions, so the
    # value shape of the 'ufl_element' (which is used to compute the
    # _value_size) will be correspond to the value size in physical
    # space. FIAT however only knows about the reference element, and
    # so the FIAT value shape of the 'element' will be the reference
    # value size. This of course only matters for elements that have
    # different physical and reference value shapes and sizes.

    return {"mappings": element.mapping(),
            "reference_value_size": _value_size(element),
            "physical_value_size": _value_size(ufl_element),
            "geometric_dimension": cell.geometric_dimension(),
            "topological_dimension": cell.topological_dimension(),
            "dofs": [L.pt_dict for L in element.dual_basis()],
            "physical_offsets": _generate_physical_offsets(ufl_element)}

def _extract_elements(element):

    new_elements = []
    if isinstance(element, (MixedElement, EnrichedElement)):
        for e in element.elements():
            new_elements += _extract_elements(e)
    else:
        new_elements.append(element)
    return new_elements

# def _num_components(element):
#     """Compute the number of components of element, like _value_size, but
#     does not support tensor elements."""
#     shape = element.value_shape()
#     if shape == ():
#         return 1
#     elif len(shape) == 1:
#         return shape[0]
#     else:
#         error("Tensor valued elements are not supported yet: %d " % shape)

def _evaluate_basis(ufl_element, element, cell):
    "Compute intermediate representation for evaluate_basis."

    # Handle Mixed and EnrichedElements by extracting 'sub' elements.
    elements = _extract_elements(element)
    offsets = _generate_reference_offsets(element) # Must check?
    mappings = element.mapping()

    # This function is evidently not implemented for TensorElements
    for e in elements:
        if len(e.value_shape()) > 1:
            return "Function not supported/implemented for TensorElements."

    # Handle QuadratureElement, not supported because the basis is only defined
    # at the dof coordinates where the value is 1, so not very interesting.
    for e in elements:
        if isinstance(e, QuadratureElement):
            return "Function not supported/implemented for QuadratureElement."

    # Initialise data with 'global' values.
    data = {"reference_value_size": _value_size(element),
            "physical_value_size": _value_size(ufl_element),
            "cellname" : cell.cellname(),
            "topological_dimension" : cell.topological_dimension(),
            "geometric_dimension" : cell.geometric_dimension(),
            "space_dimension" : element.space_dimension(),
            "needs_oriented": needs_oriented_jacobian(element),
            "max_degree": max([e.degree() for e in elements])
            }

    # Loop element and space dimensions to generate dof data.
    dof = 0
    dof_data = []
    for e in elements:
        for i in range(e.space_dimension()):
            num_components = _value_size(e)
            coefficients = []
            coeffs = e.get_coeffs()

            # Handle coefficients for vector valued basis elements
            # [Raviart-Thomas, Brezzi-Douglas-Marini (BDM)].
            if num_components > 1:
                for c in range(num_components):
                    coefficients.append(coeffs[i][c])
            else:
                coefficients.append(coeffs[i])

            dof_data.append(
              {
              "embedded_degree" : e.degree(),
              "coeffs" : coefficients,
              "num_components" : num_components,
              "dmats" : e.dmats(),
              "mapping" : mappings[dof],
              "offset" : offsets[dof],
              "num_expansion_members": e.get_num_members(e.degree())
              })

            dof += 1

    data["dof_data"] = dof_data

    return data

def _tabulate_coordinates(ufl_element, element):
    "Compute intermediate representation of tabulate_coordinates."

    if uses_integral_moments(element):
        return {}

    data = {}
    data["tdim"] = ufl_element.cell().topological_dimension()
    data["gdim"] = ufl_element.cell().geometric_dimension()
    data["points"] = [L.pt_dict.keys()[0] for L in element.dual_basis()]
    return data

def _tabulate_dofs(element, cell):
    "Compute intermediate representation of tabulate_dofs."

    if isinstance(element, SpaceOfReals):
        return None

    # Extract number of entities for each dimension for this cell
    num_entities = cellname_to_num_entities[cell.cellname()]

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
    multiple_entities =  any([sum(n > 0 for n in num_dofs) - 1
                              for num_dofs in num_dofs_per_element])
    need_offset = len(elements) > 1 or multiple_entities

    num_dofs_per_element = [e.space_dimension() for e in elements]

    # Handle global "elements"
    fakes = [isinstance(e, SpaceOfReals) for e in elements]

    return (dofs_per_element, num_dofs_per_element, num_entities, need_offset, fakes)


def _tabulate_facet_dofs(element, cell):
    "Compute intermediate representation of tabulate_facet_dofs."

    # Compute incidences
    incidence = __compute_incidence(cell.topological_dimension())

    # Get topological dimension
    D = max([pair[0][0] for pair in incidence])

    # Get the number of facets
    num_facets = cellname_to_num_entities[cell.cellname()][-2]

    # Find out which entities are incident to each facet
    incident = num_facets*[None]
    for facet in range(num_facets):
        incident[facet] = [pair[1] for pair in incidence if incidence[pair] == True and pair[0] == (D - 1, facet)]

    # Make list of dofs
    facet_dofs = []
    entity_dofs = element.entity_dofs()

    for facet in range(num_facets):
        facet_dofs += [[]]
        for dim in entity_dofs:
            for entity in entity_dofs[dim]:
                if (dim, entity) in incident[facet]:
                    facet_dofs[facet] += entity_dofs[dim][entity]
        facet_dofs[facet].sort()
    return facet_dofs

def _interpolate_vertex_values(ufl_element, element, cell):
    "Compute intermediate representation of interpolate_vertex_values."

    # Check for QuadratureElement
    for e in all_elements(element):
        if isinstance(e, QuadratureElement):
            return "Function is not supported/implemented for QuadratureElement."

    ir = {}
    ir["geometric_dimension"] = cell.geometric_dimension()
    ir["topological_dimension"] = cell.topological_dimension()

    # Check whether computing the Jacobian is necessary
    mappings = element.mapping()
    ir["needs_jacobian"] = any("piola" in m for m in mappings)
    ir["needs_oriented"] = needs_oriented_jacobian(element)

    # See note in _evaluate_dofs
    ir["reference_value_size"] = _value_size(element)
    ir["physical_value_size"] = _value_size(ufl_element)

    # Get vertices of reference cell
    cell = reference_cell(cell.cellname())
    vertices = cell.get_vertices()

    # Compute data for each constituent element
    extract = lambda values: values[values.keys()[0]].transpose()
    all_fiat_elm = all_elements(element)
    ir["element_data"] = [{
                           # See note in _evaluate_dofs
                           "reference_value_size": _value_size(e),
                           "physical_value_size": _value_size(e), # FIXME: Get from corresponding ufl element
                           "basis_values": extract(e.tabulate(0, vertices)),
                           "mapping": e.mapping()[0],
                           "space_dim": e.space_dimension()}
                          for e in all_fiat_elm]

    # FIXME: Temporary hack!
    if len(ir["element_data"]) == 1:
        ir["element_data"][0]["physical_value_size"] = ir["physical_value_size"]

    # Consistency check, related to note in _evaluate_dofs
    # This will fail for e.g. (RT1 x DG0) on a manifold
    if sum(data["physical_value_size"] for data in ir["element_data"]) != ir["physical_value_size"]:
        ir = "Failed to set physical value size correctly for subelements."
    elif sum(data["reference_value_size"] for data in ir["element_data"]) != ir["reference_value_size"]:
        ir = "Failed to set reference value size correctly for subelements."

    return ir

def _create_sub_foo(ufl_element, element_numbers):
    "Compute intermediate representation of create_sub_element/dofmap."
    return [element_numbers[e] for e in ufl_element.sub_elements()]

def _create_foo_integral(domain_type, form_data):
    "Compute intermediate representation of create_foo_integral."
    return [itg_data.domain_id for itg_data in form_data.integral_data
           if itg_data.domain_type == domain_type and isinstance(itg_data.domain_id, int)]

def _has_foo_integrals(domain_type, form_data):
    "Compute intermediate representation of has_foo_integrals."
    v = (form_data.num_sub_domains.get(domain_type,0) > 0
         or _create_default_foo_integral(domain_type, form_data) is not None)
    return bool(v)

def _create_default_foo_integral(domain_type, form_data):
    "Compute intermediate representation of create_default_foo_integral."
    itg_data = [itg_data for itg_data in form_data.integral_data
           if itg_data.domain_id == Measure.DOMAIN_ID_OTHERWISE and itg_data.domain_type == domain_type]
    ffc_assert(len(itg_data) in (0,1), "Expecting at most one default integral of each type.")
    return Measure.DOMAIN_ID_OTHERWISE if itg_data else None

#--- Utility functions ---

# FIXME: KBO: This could go somewhere else, like in UFL?
#        MSA: There is probably something related in ufl somewhere,
#        but I don't understand quite what this does.
#        In particular it does not cover sub-sub-elements? Is that a bug?
# Also look at function naming, use single '_' for utility functions.
def all_elements(element):

    if isinstance(element, MixedElement):
        return element.elements()

    return [element]

def _num_dofs_per_entity(element):
    """
    Compute list of integers representing the number of dofs
    associated with a single mesh entity.

    Example: Lagrange of degree 3 on triangle: [1, 2, 1]
    """
    entity_dofs = element.entity_dofs()
    return [len(entity_dofs[e][0]) for e in range(len(entity_dofs.keys()))]

# These two are copied from old ffc
def __compute_incidence(D):
    "Compute which entities are incident with which"

    # Compute the incident vertices for each entity
    sub_simplices = []
    for dim in range(D + 1):
        sub_simplices += [__compute_sub_simplices(D, dim)]

    # Check which entities are incident, d0 --> d1 for d0 >= d1
    incidence = {}
    for d0 in range(0, D + 1):
        for i0 in range(len(sub_simplices[d0])):
            for d1 in range(d0 + 1):
                for i1 in range(len(sub_simplices[d1])):
                    if min([v in sub_simplices[d0][i0] for v in sub_simplices[d1][i1]]) == True:
                        incidence[((d0, i0), (d1, i1))] = True
                    else:
                        incidence[((d0, i0), (d1, i1))] = False

    return incidence

def __compute_sub_simplices(D, d):
    """Compute vertices for all sub simplices of dimension d (code
    taken from Exterior)."""

    # Number of vertices
    num_vertices = D + 1

    # Special cases: d = 0 and d = D
    if d == 0:
        return [[i] for i in range(num_vertices)]
    elif d == D:
        return [range(num_vertices)]

    # Compute all permutations of num_vertices - (d + 1)
    permutations = compute_permutations(num_vertices - d - 1, num_vertices)

    # Iterate over sub simplices
    sub_simplices = []
    for i in range(len(permutations)):

        # Pick tuple i among permutations (non-incident vertices)
        remove = permutations[i]

        # Remove vertices, keeping d + 1 vertices
        vertices = [v for v in range(num_vertices) if not v in remove]
        sub_simplices += [vertices]

    return sub_simplices

def uses_integral_moments(element):
    "True if element uses integral moments for its degrees of freedom."

    integrals = set(["IntegralMoment", "FrobeniusIntegralMoment"])
    tags = set([L.get_type_tag() for L in element.dual_basis()])
    return len(integrals & tags) > 0

def needs_oriented_jacobian(element):
    # Check whether this element needs an oriented jacobian (only
    # contravariant piolas seem to need it)
    return ("contravariant piola" in element.mapping())
