# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2018.

declaration = """
extern ufc_dofmap {factory_name};
"""

factory = """
// Code for dofmap {factory_name}

{sub_dofmaps_initialization}

ufc_dofmap {factory_name} =
{{
  .signature = {signature},
  .block_size = {block_size},
  .num_global_support_dofs = {num_global_support_dofs},
  .num_element_support_dofs = {num_element_support_dofs},
  .num_sub_dofmaps = {num_sub_dofmaps},
  .sub_dofmaps = {sub_dofmaps}
}};

// End of code for dofmap {factory_name}
"""
