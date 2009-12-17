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

# Last changed: 2009-12-17

from log import debug_ir
from fiatinterface import create_fiat_element

not_implemented = None

def compute_form_ir(form, form_data, method):
    "Compute and return intermediate representation of form."

    # FIXME: Call correct method in quadrature or tensor module

    if method == "quadrature":
        return {}
    else:
        return {}

def compute_element_ir(ufl_element):
    "Compute and return intermediate representation of element."

    # Note to developers: oneliners or call a _function

    # Create FIAT element
    fiat_element = create_fiat_element(ufl_element)

    # Compute data for each function
    ir = {}
    ir["signature"] = repr(ufl_element)
    ir["cell_shape"] = ufl_element.cell().domain()
    ir["space_dimension"] = fiat_element.space_dimension()
    ir["value_rank"] = fiat_element.value_rank()
    ir["value_dimension"] = fiat_element.value_shape()
    ir["evaluate_basis"] = fiat_element.get_coeffs()
    ir["evaluate_basis_all"] = fiat_element.get_coeffs()
    ir["evaluate_basis_derivatives"] = fiat_element.get_coeffs()
    ir["evaluate_basis_derivatives_all"] = fiat_element.get_coeffs()
    ir["evaluate_dof"] = None
    ir["evaluate_dofs"] = None
    ir["interpolate_vertex_values"] = None
    ir["num_sub_elements"] = fiat_element.num_sub_elements()
    ir["create_sub_element"] = None

    debug_ir(ir, "finite_element")

    return ir

def compute_dofmap_ir(ufl_element):
    "Compute and return intermediate representation of dofmap."

    # Note to developers: oneliners or call a _function

    # Create FIAT element
    fiat_element = create_fiat_element(ufl_element)

    # Compute data for each function
    ir = {}
    ir["dof_map* create_sub_dof_map"] = not_implemented
    ir["init_cell"] = not_implemented
    ir["init_cell_finalize"] = not_implemented
    ir["init_mesh"] = not_implemented
    ir["local_dimension"] = fiat_element.space_dimension()
    ir["geometric_dimension"] = fiat_element.geometric_dimension()
    ir["global_dimension"] = not_implemented
    ir["max_local_dimension"] = fiat_element.space_dimension()
    ir["needs_mesh_entities"] = _needs_mesh_entities(fiat_element)
    ir["num_entity_dofs"] = _num_dofs_per_dim(fiat_element)
    ir["num_facet_dofs"] = _num_facet_dofs(fiat_element)
    ir["num_sub_dof_maps"] =  fiat_element.num_sub_elements()
    ir["signature"] = "FFC dofmap for " + repr(ufl_element)
    ir["tabulate_dofs"] = not_implemented
    ir["tabulate_facet_dofs"] = not_implemented
    ir["tabulate_entity_dofs"] = not_implemented
    ir["tabulate_coordinates"] = not_implemented

    debug_ir(ir, "dofmap")

    return ir

#--- Utility functions ---

def _needs_mesh_entities(fiat_element):
    """
    Returns tuple (bool, ..., bool) with item i being true/false
    depending on whether mesh entity of dimension i is needed.
    """

    return [d > 0 for d in _num_dofs_per_dim(fiat_element)]

def _num_dofs_per_dim(element):
    """Compute the number of dofs associated with each topological
    dimension.  Currently only handles non-mixed elements.
    """

    return [sum(len(dof_indices) for dof_indices in entity_dof.values())
            for (i, entity_dof) in element.entity_dofs().iteritems()]


def _num_facet_dofs(element):
    "Compute the number of dofs on associated with each cell facet."

    dim = element.geometric_dimension()
    num_facet_entities = {1: (1, 0), 2: (2, 1, 0), 3: (3, 3, 1, 0)}[dim]
    entity_dofs = element.entity_dofs()

    return sum(len(entity_dofs[entity][0])*num
               for (entity, num) in enumerate(num_facet_entities))
