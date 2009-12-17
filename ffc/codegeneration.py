"""
This module implements the generation of C++ code for the body of each
UFC function from an (optimized) intermediate representation (OIR).
"""

__author__ = "Anders Logg (logg@simula.no) and friends"
__date__ = "2009-12-16"
__copyright__ = "Copyright (C) 2009 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2009-12-17

from log import debug_code
from cpp import format

def generate_element_code(ir):
    "Generate code for finite element from intermediate representation."

    # Prefetch formatting to speedup code generation
    ret = format["return"]

    # Generate code
    code = {}
    code["classname"] = "FooClass"
    code["members"] = ""
    code["constructor"] = ""
    code["destructor"] = ""
    code["signature"] = ret('"%s"' % ir["signature"])
    code["cell_shape"] = ""
    code["space_dimension"] = ret(ir["space_dimension"])
    code["value_rank"] = ret(ir["value_rank"])
    code["value_dimension"] = ""
    code["evaluate_basis"] = ""
    code["evaluate_basis_all"] = ""
    code["evaluate_basis_derivatives"] = ""
    code["evaluate_basis_derivatives_all"] = ""
    code["evaluate_dof"] = ""
    code["evaluate_dofs"] = ""
    code["interpolate_vertex_values"] = ""
    code["num_sub_elements"] = ret(ir["num_sub_elements"])
    code["create_sub_element"] = ""

    debug_code(code, "finite_element")

    return code

def generate_dofmap_code(ir):
    "Generate code for dofmap from intermediate representation."

    # Prefetch formatting to speedup code generation
    ret = format["return"]

    # Generate code
    code = {}
    code["classname"] = "FooClass"
    code["members"] = ""
    code["constructor"] = ""
    code["destructor"] = ""
    code["signature"] = ret(ir["signature"])
    code["needs_mesh_entities"] = ""
    code["init_mesh"] = ""
    code["init_cell"] = ""
    code["init_cell_finalize"] = ""
    code["global_dimension"] = ret("_global_dimension")
    code["local_dimension"] = ret(ir["local_dimension"])
    code["max_local_dimension"] = ret(ir["max_local_dimension"])
    code["geometric_dimension"] = ret(ir["geometric_dimension"])
    code["num_facet_dofs"] = ret(ir["local_dimension"])
    code["num_entity_dofs"] = ""
    code["tabulate_dofs"] = ""
    code["tabulate_facet_dofs"] = ""
    code["tabulate_entity_dofs"] = ""
    code["tabulate_coordinates"] = ""
    code["num_sub_dof_maps"] = ret(ir["num_sub_dof_maps"])
    code["create_sub_dof_map"] = ""

    debug_code(code, "dofmap")

    return code
