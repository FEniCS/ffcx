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

    code = {}
    code["signature"] = format["return"]('"%s"' % ir["signature"])
    code["cell_shape"] = ""
    code["space_dimension"] = format["return"](ir["space_dimension"])
    code["value_rank"] = format["return"](ir["value_rank"])
    code["value_dimension"] = ""
    code["evaluate_basis"] = ""
    code["evaluate_basis_all"] = ""
    code["evaluate_basis_derivatives"] = ""
    code["evaluate_basis_derivatives_all"] = ""
    code["evaluate_dof"] = ""
    code["evaluate_dofs"] = ""
    code["interpolate_vertex_values"] = ""
    code["num_sub_elements"] = ""
    code["create_sub_element"] = ""

    debug_code(code, "finite_element")

    return code

def generate_dofmap_code(ir):
    "Generate code for dofmap from intermediate representation."

    print ir

    code = {}
    code["signature"] = format["return"](ir["signature"])
    code["needs_mesh_entities"] = ""
    code["init_mesh"] = ""
    code["init_cell"] = ""
    code["init_cell_finalize"] = ""
    code["global_dimension"] = ""
    code["local_dimension"] = format["return"](ir["local_dimension"])
    code["max_local_dimension"] = ""
    code["geometric_dimension"] = format["return"](ir["geometric_dimension"])
    code["num_facet_dofs"] = format["return"](ir["local_dimension"])
    code["num_entity_dofs"] = ""
    code["tabulate_dofs"] = ""
    code["tabulate_facet_dofs"] = ""
    code["tabulate_entity_dofs"] = ""
    code["tabulate_coordinates"] = ""
    code["num_sub_dof_maps"] = ""
    code["dof_map* create_sub_dof_map"] = ""

    debug_code(code, "dofmap")

    return code
