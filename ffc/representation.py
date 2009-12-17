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
    ir["signature"] = "FFC dofmap for " + repr(ufl_element)
    ir["needs_mesh_entities"] = None
    ir["init_mesh"] = None
    ir["init_cell"] = None
    ir["init_cell_finalize"] = None
    ir["global_dimension"] = None
    ir["local_dimension"] = fiat_element.space_dimension()
    ir["max_local_dimension"] = fiat_element.space_dimension()
    ir["geometric_dimension"] = fiat_element.geometric_dimension()
    ir["num_facet_dofs"] = _num_facet_dofs(fiat_element)
    ir["num_entity_dofs"] = None
    ir["tabulate_dofs"] = None
    ir["tabulate_facet_dofs"] = None
    ir["tabulate_entity_dofs"] = None
    ir["tabulate_coordinates"] = None
    ir["num_sub_dof_maps"] = fiat_element.num_sub_elements()
    ir["dof_map* create_sub_dof_map"] = None

    debug_ir(ir, "dofmap")

    return {}

#--- Utility functions ---

def _num_facet_dofs(fiat_element):
    "Compute the number of dofs on each cell facet."

    # FIXME: Seems to need updating for new FIAT interface
    return 1

    num_facet_entities = {"interval": (1, 0),
                          "triangle": (2, 1, 0),
                          "tetrahedron": (3, 3, 1, 0)}
    num_dofs_per_dim = _num_dofs_per_dim(fiat_element)
    num_facet_dofs = 0
    for dim in range(len(num_dofs_per_dim)):
        num_facet_dofs += num_facet_entities[cell_shape][dim]*num_dofs_per_dim[dim]

    return num_facet_dofs

def _num_dofs_per_dim(fiat_element):
    "Compute the number of dofs associated with each topological dimension."
    num_dofs_per_dim = []
    for sub_entity_dofs in fiat_element.entity_dofs():
        sub_num_dofs_per_dim = {}
        for dim in sub_entity_dofs:
            num_dofs = [len(sub_entity_dofs[dim][e]) for e in sub_entity_dofs[dim]]
            if dim in sub_num_dofs_per_dim:
                sub_num_dofs_per_dim[dim] += pick_first(num_dofs)
            else:
                sub_num_dofs_per_dim[dim] = pick_first(num_dofs)
        num_dofs_per_dim += [sub_num_dofs_per_dim]
    return num_dofs_per_dim
