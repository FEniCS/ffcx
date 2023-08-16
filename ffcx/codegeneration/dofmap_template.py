# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2018.

declaration = """
extern ufcx_dofmap {factory_name};
"""

factory = """
// Code for dofmap {factory_name}

{sub_dofmaps_initialization}

{entity_dofs_init}

{entity_dof_offsets_init}

{entity_closure_dofs_init}

{entity_closure_dof_offsets_init}

ufcx_dofmap {factory_name} =
{{
  .signature = {signature},
  .num_global_support_dofs = {num_global_support_dofs},
  .num_element_support_dofs = {num_element_support_dofs},
  .block_size = {block_size},
  .entity_dofs = {entity_dofs},
  .entity_dof_offsets = {entity_dof_offsets},
  .entity_closure_dofs = {entity_closure_dofs},
  .entity_closure_dof_offsets = {entity_closure_dof_offsets},
  .num_sub_dofmaps = {num_sub_dofmaps},
  .sub_dofmaps = {sub_dofmaps}
}};

// End of code for dofmap {factory_name}
"""
