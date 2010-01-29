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

__author__ = "Anders Logg (logg@simula.no) and friends"
__date__ = "2009-12-16"
__copyright__ = "Copyright (C) 2009-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Marie E. Rognes 2010
# Modified by Kristian B. Oelgaard 2010
# Last changed: 2010-01-29

# Python modules
from itertools import chain

# Import UFL
import ufl

# FFC modules
from ffc.utils import compute_permutations, product
from ffc.log import info, error, begin, end, debug_ir, ffc_assert, warning
from ffc.fiatinterface import create_element, entities_per_dim, reference_cell
from ffc.mixedelement import MixedElement
from ffc.quadratureelement import QuadratureElement

# FFC specialized representation modules
from ffc import quadrature
from ffc import tensor

not_implemented = None

def compute_ir(analysis, options):
    "Compute intermediate representation."

    begin("Compiler stage 2: Computing intermediate representation")

    # Extract data from analysis
    form_and_data, elements, element_map = analysis

    # Compute representation of elements
    info("Computing representation of %d elements" % len(elements))
    ir_elements = [_compute_element_ir(e, i, element_map) for (i, e) in enumerate(elements)]

    # Compute representation of dofmaps
    info("Computing representation of %d dofmaps" % len(elements))
    ir_dofmaps = [_compute_dofmap_ir(e, i, element_map) for (i, e) in enumerate(elements)]

    # Compute and flatten representation of integrals
    info("Computing representation of integrals")
    irs = [_compute_integral_ir(f, d, i) for (i, (f, d)) in enumerate(form_and_data)]
    ir_integrals = [ir for ir in chain(*irs) if not ir is None]

    # Compute representation of forms
    info("Computing representation of forms")
    ir_forms = [_compute_form_ir(f, d, i, element_map) for (i, (f, d)) in enumerate(form_and_data)]

    end()

    return ir_elements, ir_dofmaps, ir_integrals, ir_forms

def _compute_element_ir(ufl_element, element_id, element_map):
    "Compute intermediate representation of element."

    # Create FIAT element
    element = create_element(ufl_element)
    cell = ufl_element.cell()

    # Store id
    ir = {"id": element_id}

    # Compute data for each function
    ir["signature"] = repr(ufl_element)
    ir["cell_shape"] = cell.domain()
    ir["space_dimension"] = element.space_dimension()
    ir["value_rank"] = len(ufl_element.value_shape())
    ir["value_dimension"] = ufl_element.value_shape()
    ir["evaluate_basis"] = _verify_evaluate_basis(_evaluate_basis(ufl_element, element))
    ir["evaluate_dof"] = _evaluate_dof(element, cell)
    ir["interpolate_vertex_values"] = _interpolate_vertex_values(element, cell)
    ir["num_sub_elements"] = ufl_element.num_sub_elements()
    ir["create_sub_element"] = _create_sub_foo(ufl_element, element_map)

    #debug_ir(ir, "finite_element")

    return ir

def _compute_dofmap_ir(ufl_element, element_id, element_map):
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
    ir["needs_mesh_entities"] = [d > 0 for d in num_dofs_per_entity]
    ir["init_mesh"] = num_dofs_per_entity
    ir["init_cell"] = None
    ir["init_cell_finalize"] = None
    ir["global_dimension"] = None
    ir["local_dimension"] = element.space_dimension()
    ir["max_local_dimension"] = element.space_dimension()
    ir["geometric_dimension"] = cell.geometric_dimension()
    ir["num_facet_dofs"] = len(facet_dofs[0])
    ir["num_entity_dofs"] = num_dofs_per_entity
    ir["tabulate_dofs"] = _tabulate_dofs(element, cell)
    ir["tabulate_facet_dofs"] = facet_dofs
    ir["tabulate_entity_dofs"] = (element.entity_dofs(), num_dofs_per_entity)
    ir["tabulate_coordinates"] = _tabulate_coordinates(element)
    ir["num_sub_dof_maps"] = ufl_element.num_sub_elements()
    ir["create_sub_dof_map"] = _create_sub_foo(ufl_element, element_map)

    #debug_ir(ir, "dofmap")

    return ir

def _compute_integral_ir(form, form_data, form_id):
    "Compute intermediate represention of form integrals."

    irs = []

    # Iterate over integrals
    for (domain_type, domain_id, integrals, metadata) in form_data.integral_data:

        # Select representation
        if metadata["representation"] == "quadrature":
            r = quadrature
        elif metadata["representation"] == "tensor":
            r = tensor
        else:
            error("Unknown representation: " + str(metadata["representation"]))

        # Compute representation
        ir = r.compute_integral_ir(domain_type,
                                   domain_id,
                                   integrals,
                                   metadata,
                                   form_data,
                                   form_id)

        # Append representation
        irs.append(ir)

    return irs

def _compute_form_ir(form, form_data, form_id, element_map):
    "Compute intermediate representation of form."

    # Store id
    ir = {"id": form_id}

    # Compute common data
    ir["classname"] = "FooForm"
    ir["members"] = not_implemented
    ir["constructor"] = not_implemented
    ir["destructor"] = not_implemented
    ir["signature"] = repr(form)
    ir["rank"] = form_data.rank
    ir["num_coefficients"] = form_data.num_coefficients
    ir["num_cell_integrals"] = form_data.num_cell_domains
    ir["num_exterior_facet_integrals"] = form_data.num_exterior_facet_domains
    ir["num_interior_facet_integrals"] = form_data.num_interior_facet_domains
    ir["create_finite_element"] = [element_map[e] for e in form_data.elements]
    ir["create_dof_map"] = [element_map[e] for e in form_data.elements]
    ir["create_cell_integral"] = _create_foo_integral("cell", form_data)
    ir["create_exterior_facet_integral"] = _create_foo_integral("exterior_facet", form_data)
    ir["create_interior_facet_integral"] = _create_foo_integral("interior_facet", form_data)

    return ir

#--- Computation of intermediate representation for non-trivial functions ---

def _evaluate_basis(ufl_element, fiat_element):
    "Compute intermediate representation for evaluate_basis."

    # Handle MixedElements recursively.
    # TODO: KBO: Is this OK if done consistently in FFC (for TensorElements)?
    if isinstance(ufl_element, ufl.MixedElement):
        data = []
        for element in ufl_element.sub_elements():
            data += _evaluate_basis(element, create_element(element))
        return data

    # Handle QuadratureElement, not supported because the basis is only defined
    # at the dof coordinates where the value is 1, so not very interesting.
    if ufl_element.family() == "Quadrature":
        return [not_implemented]

    # TODO: KBO: Remove if never triggered.
    ffc_assert(fiat_element.get_nodal_basis().get_embedded_degree() == \
               ufl_element.degree(),\
               "Degrees do not match: %s, %s" % (repr(fiat_element), repr(ufl_element)))

    # Get mapping for element, must be the same for all DOFs for evaluate_basis
    # to work in its current implementation
    mapping = fiat_element.mapping()[0]
    ffc_assert(all(mapping == m for m in fiat_element.mapping()),\
               "Mapping is not the same for all dofs of this element: %s" % str(fiat_element))
    data = {
          "value_shape" : fiat_element.value_shape(),
          "embedded_degree" : ufl_element.degree(),
          "cell_domain" : ufl_element.cell().domain(),
          "coeffs" : fiat_element.get_coeffs(),
          "family" : ufl_element.family(),
          "mapping" : mapping,
          "space_dimension" : fiat_element.space_dimension(),
          "topological_dimension" : ufl_element.cell().topological_dimension(),
          "geometric_dimension" : ufl_element.cell().geometric_dimension(),
          "dmats" : fiat_element.get_nodal_basis().get_dmats()
          }

    data["num_expansion_members"] = \
      fiat_element.get_nodal_basis().get_expansion_set().get_num_members(data["embedded_degree"])

    return [data]

def _verify_evaluate_basis(data_list):
    "Some safety checks."

    # If we have a QuadratureElement, or if there is one present in the MixedElement
    # we cannot create evaluate_basis
    for data in data_list:
        if data is not_implemented:
            return "Function not supported/implemented for QuadratureElement."


    # FIXME: KBO: If UFL makes sure that the below is always true, we can delete this function.
    # Get the element cell domain and check if it is the same for all elements.
    element_cell_domain = data_list[0]["cell_domain"]
    ffc_assert(all(element_cell_domain == data["cell_domain"] for data in data_list),\
               "The element cell domain must be the same for all sub elements: " + repr(data_list))

    # Get the element geometric dimension and check if it is the same for all elements.
    geometric_dimension = data_list[0]["geometric_dimension"]
    ffc_assert(all(geometric_dimension == data["geometric_dimension"] for data in data_list),\
               "The geometric dimension must be the same for all sub elements: " + repr(data_list))

    # Get the element topological dimension and check if it is the same for all elements.
    topological_dimension = data_list[0]["topological_dimension"]
    ffc_assert(all(topological_dimension == data["topological_dimension"] for data in data_list),\
               "The topological dimension must be the same for all sub elements: " + repr(data_list))
    return data_list

# FIXME: Move to FiniteElement/MixedElement
def _value_size(element):
    "Compute value size of element."
    shape = element.value_shape()
    if shape == ():
        return 1
    else:
        return product(shape)

def _evaluate_dof(element, cell):
    "Compute intermediate representation of evaluate_dof."

    # Generate offsets: i.e value offset for each basis function
    offsets = []
    offset = 0
    for e in all_elements(element):
        offsets += [offset]*e.space_dimension()
        offset += _value_size(e) # AL: change here, is that correct?
    return {"mappings": element.mapping(),
            "value_size": _value_size(element),
            "cell_dimension": cell.geometric_dimension(),
            "dofs": [L.pt_dict for L in element.dual_basis()],
            "offsets": offsets}

def _tabulate_coordinates(element):
    "Compute intermediate representation of tabulate_coordinates."
    if uses_integral_moments(element):
        return None
    return [L.pt_dict.keys()[0] for L in element.dual_basis()]

def _tabulate_dofs(element, cell):
    "Compute intermediate representation of tabulate_dofs."

    # Extract number of enties for each dimension for this cell
    num_entities = entities_per_dim[cell.geometric_dimension()]

    # Extract number of dofs per entity for each element
    elements = all_elements(element)
    num_dofs_per_element = [_num_dofs_per_entity(e) for e in elements]

    # Check whether we need offset
    multiple_entities =  any([sum(n > 0 for n in num_dofs) - 1
                              for num_dofs in num_dofs_per_element])
    need_offset = len(elements) > 1 or multiple_entities

    return (num_dofs_per_element, num_entities, need_offset)


def _tabulate_facet_dofs(element, cell):
    "Compute intermediate representation of tabulate_facet_dofs."

    # Compute incidences
    incidence = __compute_incidence(cell.geometric_dimension())

    # Get topological dimension
    D = max([pair[0][0] for pair in incidence])

    # Get the number of facets
    num_facets = cell.num_facets()

    # Find out which entities are incident to each facet
    incident = num_facets*[[]]
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

def _interpolate_vertex_values(element, cell):
    "Compute intermediate representation of interpolate_vertex_values."

    # Check for QuadratureElement
    for e in all_elements(element):
        if isinstance(e, QuadratureElement):
            return "Function is not supported/implemented for QuadratureElement."

    ir = {}
    ir["cell_dim"] = cell.geometric_dimension()

    # Check whether computing the Jacobian is necessary
    mappings = element.mapping()
    ir["needs_jacobian"] = any("piola" in m for m in mappings)

    # Get vertices of reference cell
    cell = reference_cell(cell.domain())
    vertices = cell.get_vertices()

    # Compute data for each constituent element
    extract = lambda values: values[values.keys()[0]].transpose()
    ir["element_data"] = [{"value_size": _value_size(e),
                           "basis_values": extract(e.tabulate(0, vertices)),
                           "mapping": e.mapping()[0],
                           "space_dim": e.space_dimension()}
                          for e in all_elements(element)]
    return ir

def _create_sub_foo(ufl_element, element_map):
    "Compute intermediate representation of create_sub_element/dofmap."
    return [element_map[e] for e in ufl_element.sub_elements()]

def _create_foo_integral(domain_type, form_data):
    "Compute intermediate representation of create_foo_integral."
    return [domain_id for (_domain_type, domain_id, integrals, metadata) in
           form_data.integral_data if _domain_type == domain_type]

#--- Utility functions ---

# FIXME: KBO: This could go somewhere else, like in UFL?
# Also look at function naming, use single '_' for utility functions.
def all_elements(element):
    try:
        return element.elements()
    except:
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
    "Compute vertices for all sub simplices of dimension d (code taken from Exterior)"

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

    integrals = set(["IntegralMoment", "FrobeniusIntegralMoment"])
    tags = set([L.get_type_tag() for L in element.dual_basis()])
    return len(integrals & tags) > 0
