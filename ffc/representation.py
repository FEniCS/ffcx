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
__copyright__ = "Copyright (C) 2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2009-12-16

def form_representation(form, form_data, method):
    "Compute and return intermediate representation of form."

    # FIXME: Call correct method in quadrature or tensor module

    if method == "quadrature":
        return {}
    else:
        return {}

def element_representation(ufl_element):
    "Compute and return intermediate representation of element."

    ir = {}

    ir["signature"] = repr(ufl_element)
    ir["cell_shape"] = None
    ir["space_dimension"] = None
    ir["value_rank"] = None
    ir["value_dimension"] = None
    ir["evaluate_basis"] = None
    ir["evaluate_basis_all"] = None
    ir["evaluate_basis_derivatives"] = None
    ir["evaluate_basis_derivatives_all"] = None
    ir["evaluate_dof"] = None
    ir["evaluate_dofs"] = None
    ir["interpolate_vertex_values"] = None
    ir["num_sub_elements"] = None
    ir["create_sub_element"] = None

    return ir

def dofmap_representation(ufl_element):
    "Compute and return intermediate representation of dofmap."

    ir = {}

    ir["signature"] = None
    ir["needs_mesh_entities"] = None
    ir["init_mesh"] = None
    ir["init_cell"] = None
    ir["init_cell_finalize"] = None
    ir["global_dimension"] = None
    ir["local_dimension"] = None
    ir["max_local_dimension"] = None
    ir["geometric_dimension"] = None
    ir["num_facet_dofs"] = None
    ir["num_entity_dofs"] = None
    ir["tabulate_dofs"] = None
    ir["tabulate_facet_dofs"] = None
    ir["tabulate_entity_dofs"] = None
    ir["tabulate_coordinates"] = None
    ir["num_sub_dof_maps"] = None
    ir["dof_map* create_sub_dof_map"] = None

    return {}
