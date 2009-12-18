"""
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

# Last changed: 2009-12-18

from ufl.finiteelement import FiniteElement as UFLFiniteElement
from ufl.algorithms.analysis import extract_sub_elements

from log import debug_ir
from fiatinterface import create_element
from mixedelement import MixedElement

not_implemented = None

def compute_form_ir(form, form_data, method):
    "Compute and return intermediate representation of form."

    # FIXME: Call correct method in quadrature or tensor module

    if method == "quadrature":
        return {}
    else:
        return {}

def compute_elements_ir(ufl_element):
    """
    Compute and return intermediate representation for every element
    in the 'element tree' spanned by the ufl_element.
    """
    ir = {}
    elements = extract_sub_elements(ufl_element)
    for element in elements:
        ir[element] = compute_element_ir(ufl_element)
    return ir

def compute_element_ir(ufl_element):
    "Compute and return intermediate representation of element."

    # Note to developers: oneliners or call a _function

    # Create FIAT element
    element = create_element(ufl_element)

    ir = {}
    # Compute data for each function
    ir["signature"] = repr(ufl_element)
    ir["cell_shape"] = ufl_element.cell().domain()
    ir["space_dimension"] = element.space_dimension()
    ir["value_rank"] = len(ufl_element.value_shape())
    ir["value_dimension"] = ufl_element.value_shape()
    ir["evaluate_basis"] = element
    ir["evaluate_basis_all"] = not_implemented #element.get_coeffs()
    ir["evaluate_basis_derivatives"] = element
    ir["evaluate_basis_derivatives_all"] = not_implemented #element.get_coeffs()
    ir["evaluate_dof"] = None
    ir["evaluate_dofs"] = None
    ir["interpolate_vertex_values"] = None
    ir["num_sub_elements"] = _num_sub_elements(ufl_element)
    ir["create_sub_element"] = None

    debug_ir(ir, "finite_element")

    return ir

def compute_dofmap_ir(ufl_element):
    "Compute and return intermediate representation of dofmap."

    # Note to developers: oneliners or call a _function

    # Create FIAT element
    element = create_element(ufl_element)

    # Precompute frequently used list: number of dofs per mesh entity:
    num_dofs_per_entity = _num_dofs_per_entity(element)

    # Compute data for each function
    ir = {}
    ir["dof_map* create_sub_dof_map"] = not_implemented
    ir["init_cell"] = None
    ir["init_cell_finalize"] = None
    ir["init_mesh"] = num_dofs_per_entity
    ir["local_dimension"] = element.space_dimension()
    ir["geometric_dimension"] = ufl_element.cell().domain()
    ir["global_dimension"] = None
    ir["max_local_dimension"] = element.space_dimension()
    ir["needs_mesh_entities"] = [d > 0 for d in num_dofs_per_entity]
    ir["num_entity_dofs"] = _num_dofs_per_dim(element)
    ir["num_facet_dofs"] = _num_facet_dofs(element)
    ir["num_sub_dof_maps"] =  _num_sub_elements(ufl_element)
    ir["signature"] = "FFC dofmap for " + repr(ufl_element)
    ir["tabulate_dofs"] = not_implemented
    ir["tabulate_facet_dofs"] = not_implemented
    ir["tabulate_entity_dofs"] = not_implemented
    ir["tabulate_coordinates"] = not_implemented

    debug_ir(ir, "dofmap")

    return ir

#--- Utility functions ---

def _num_sub_elements(ufl_element):
    """
    Return the number of sub_elements for a ufl_element. Would work
    better if an UFL element had a member num_sub_elements()
    """
    if isinstance(ufl_element, UFLFiniteElement):
        return 0
    else:
        return len(ufl_element.sub_elements())

def _num_dofs_per_dim(element):
    """
    Compute the number of dofs associated with each topological
    dimension.  Currently only handles non-mixed elements.

    Example: Lagrange of degree 3 on triangle:  [3, 6, 1]
    """

    entity_dofs = element.entity_dofs()
    return [sum(len(dof_indices) for dof_indices in entity_dofs[e].values())
            for e in range(len(entity_dofs.keys()))]

def _num_dofs_per_entity(element):
    """
    Compute list of integers representing the number of dofs
    associated with a single mesh entity.

    Example: Lagrange of degree 3 on triangle: [1, 2, 1]
    """
    entity_dofs = element.entity_dofs()
    return [len(entity_dofs[e][0]) for e in range(len(entity_dofs.keys()))]

def _num_facet_dofs(element):
    """
    Compute the number of dofs on associated with each cell facet.

    Example: Lagrange of degree 3 on triangle: 4
    """

    entity_dofs = element.entity_dofs()
    dim = len(entity_dofs.keys()) - 1
    num_facet_entities = {1: (1, 0), 2: (2, 1, 0), 3: (3, 3, 1, 0)}[dim]

    return sum(len(entity_dofs[entity][0])*num
               for (entity, num) in enumerate(num_facet_entities))
