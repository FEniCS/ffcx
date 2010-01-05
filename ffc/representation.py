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
__copyright__ = "Copyright (C) 2009 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Marie E. Rognes 2010
# Modified by Kristian B. Oelgaard 2010
# Last changed: 2010-01-04

# UFL modules
from ufl.finiteelement import FiniteElement as UFLFiniteElement
from ufl.finiteelement import MixedElement as UFLMixedElement

# FFC modules
from ffc.utils import compute_permutations
from ffc.log import info, error, begin, end, debug_ir, ffc_assert
from ffc.fiatinterface import create_element, entities_per_dim
from ffc.mixedelement import MixedElement


# FFC specialized representation modules
#from ffc.quadrature import QuadratureRepresentation
from ffc.tensor import TensorRepresentation

not_implemented = None

def compute_ir(form, form_data, options):
    """
    Compute and return intermediate representation of form, including
    all associated elements, dofmaps and integrals.
    """

    begin("Compiler stage 2: Computing intermediate representation")

    # Compute representation of elements, dofmaps, forms and integrals
    ir_elements  = [compute_element_ir(e) for e in form_data.unique_sub_elements]
    ir_dofmaps   = [compute_dofmap_ir(e) for e in form_data.unique_sub_elements]
    ir_integrals = compute_integrals_ir(form, form_data, options)
    ir_form      = compute_form_ir(form, form_data)

    end()

    return ir_form, ir_elements, ir_dofmaps, ir_integrals

def compute_element_ir(ufl_element):
    "Compute and return intermediate representation of element."

    info("Computing element representation")

    # Create FIAT element
    element = create_element(ufl_element)

    # Compute data for each function
    ir = {}
    ir["signature"] = repr(ufl_element)
    ir["cell_shape"] = ufl_element.cell().domain()
    ir["space_dimension"] = element.space_dimension()
    ir["value_rank"] = len(ufl_element.value_shape())
    ir["value_dimension"] = ufl_element.value_shape()
    ir["evaluate_basis"] = _evaluate_basis_data(ufl_element, element)
    ir["evaluate_basis_all"] = not_implemented #element.get_coeffs()
    ir["evaluate_basis_derivatives"] = element
    ir["evaluate_basis_derivatives_all"] = not_implemented #element.get_coeffs()
    ir["evaluate_dof"] = _evaluate_dof_ir(element, ufl_element.cell())
    ir["evaluate_dofs"] = None
    ir["interpolate_vertex_values"] = None
    ir["num_sub_elements"] = _num_sub_elements(ufl_element)
    ir["create_sub_element"] = None

    debug_ir(ir, "finite_element")

    return ir

def compute_dofmap_ir(ufl_element):
    "Compute and return intermediate representation of dofmap."

    info("Computing dofmap representation")

    # Create FIAT element
    element = create_element(ufl_element)
    cell = ufl_element.cell()

    # Precompute repeatedly used items
    num_dofs_per_entity = _num_dofs_per_entity(element)
    facet_dofs = _tabulate_facet_dofs_ir(element, cell)

    # Get list of subelements
    if _num_sub_elements(ufl_element) == 1:
        sub_elements = [element]
    else:
        sub_elements = element._elements

    # Compute data for each function
    ir = {}
    ir["dof_map* create_sub_dof_map"] = not_implemented
    ir["init_cell"] = None
    ir["init_cell_finalize"] = None
    ir["init_mesh"] = num_dofs_per_entity
    ir["local_dimension"] = element.space_dimension()
    ir["geometric_dimension"] = cell.geometric_dimension()
    ir["global_dimension"] = None
    ir["max_local_dimension"] = element.space_dimension()
    ir["needs_mesh_entities"] = [d > 0 for d in num_dofs_per_entity]
    ir["num_entity_dofs"] = num_dofs_per_entity
    ir["num_facet_dofs"] = len(facet_dofs[0])
    ir["num_sub_dof_maps"] =  _num_sub_elements(ufl_element)
    ir["signature"] = "FFC dofmap for " + repr(ufl_element)
    ir["tabulate_dofs"] = _tabulate_dofs_ir(sub_elements, cell)
    ir["tabulate_facet_dofs"] = facet_dofs
    ir["tabulate_entity_dofs"] = not_implemented
    ir["tabulate_coordinates"] = not_implemented

    debug_ir(ir, "dofmap")

    return ir

def compute_integrals_ir(form, form_data, options):
    "Compute and return intermediate represention of integrals."

    # FIXME: Handle multiple representations here
    ir = TensorRepresentation(form, form_data)

    return ir

def compute_form_ir(form, form_data):
    "Compute and return intermediate representation of form."

    info("Computing form representation")

    # Compute common data
    ir = {}
    ir["classname"] = "FooForm"
    ir["members"] = not_implemented
    ir["constructor"] = not_implemented
    ir["destructor"] = not_implemented
    ir["signature"] = repr(form)
    ir["rank"] = form_data.rank
    ir["num_coefficients"] = form_data.num_coefficients
    ir["num_cell_integrals"] = form_data.num_cell_integrals
    ir["num_exterior_facet_integrals"] = form_data.num_exterior_facet_integrals
    ir["num_interior_facet_integrals"] = form_data.num_interior_facet_integrals
    ir["create_finite_element"] = not_implemented
    ir["create_dof_map"] = not_implemented
    ir["create_cell_integral"] = not_implemented
    ir["create_exterior_facet_integral"] = not_implemented
    ir["create_interior_facet_integral"] = not_implemented

    return ir

#--- Utility functions -

def _num_sub_elements(ufl_element):
    """
    Return the number of sub_elements for a ufl_element. Would work
    better if an UFL element had a member num_sub_elements()
    """
    if isinstance(ufl_element, UFLFiniteElement):
        return 1
    else:
        return len(ufl_element.sub_elements())


def _num_dofs_per_entity(element):
    """
    Compute list of integers representing the number of dofs
    associated with a single mesh entity.

    Example: Lagrange of degree 3 on triangle: [1, 2, 1]
    """
    entity_dofs = element.entity_dofs()
    return [len(entity_dofs[e][0]) for e in range(len(entity_dofs.keys()))]


def _evaluate_dof_ir(element, cell):

    if element.value_shape() == ():
        value_dim = 1
    else:
        value_dim = element.value_shape()[0]

    # Get list of subelements
    if isinstance(element, MixedElement):
        sub_elements = element._elements
    else:
        sub_elements = [element]

    # Generate offsets
    offsets = []
    a = 0
    for e in sub_elements:
        if e.value_shape() == ():
            sub_value_dim = 1
        else:
            sub_value_dim = e.value_shape()[0]
        offsets += [a]*e.space_dimension()
        a += sub_value_dim

    return {"mappings": list(element.mapping()),
            "value_dim": value_dim,
            "cell_dimension": cell.geometric_dimension(),
            "dofs": [L.pt_dict for L in element.dual_basis()],
            "offsets": offsets}


def _tabulate_dofs_ir(sub_elements, cell):
    ""

    return [{"entites_per_dim": entities_per_dim[cell.geometric_dimension()],
             "num_dofs_per_entity": _num_dofs_per_entity(e)}
            for e in sub_elements]


def _tabulate_facet_dofs_ir(element, cell):
    ""

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

    return facet_dofs


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

def _evaluate_basis_data(ufl_element, fiat_element):
    "Helper function to extract relevant data for evaluate_basis* functions."

    # FIXME: KBO: No support for Mixed-, Vector-, TensorElement
    if isinstance(ufl_element, UFLMixedElement):
        return not_implemented


    # TODO: KBO: Remove if never triggered.
    ffc_assert(fiat_element.get_nodal_basis().get_embedded_degree() == \
               ufl_element.degree(),\
               "Degrees do not match: %s, %s" % (repr(fiat_element), repr(ufl_element)))
    # FIXME: KBO: Does not support mixed elements yet.
    data = {
          "num_sub_elements" : _num_sub_elements(ufl_element),
          "value_shape" : fiat_element.value_shape(),
          "embedded_degree" : ufl_element.degree(),
          "cell_domain" : ufl_element.cell().domain(),
          "coeffs" : fiat_element.get_coeffs(),
          "value_rank" : len(ufl_element.value_shape())
          }

    data["num_expansion_members"] = \
      fiat_element.get_nodal_basis().get_expansion_set().get_num_members(data["embedded_degree"])

    return data

