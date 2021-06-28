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

void tabulate_entity_dofs_{factory_name}(int* restrict dofs, int d, int i)
{{
{tabulate_entity_dofs}
}}

void tabulate_entity_closure_dofs_{factory_name}(int* restrict dofs, int d, int i)
{{
{tabulate_entity_closure_dofs}
}}

ufc_dofmap {factory_name} =
{{
  .signature = {signature},
  .num_global_support_dofs = {num_global_support_dofs},
  .num_element_support_dofs = {num_element_support_dofs},
  .block_size = {block_size},
  .num_entity_dofs[0] = {num_entity_dofs[0]},
  .num_entity_dofs[1] = {num_entity_dofs[1]},
  .num_entity_dofs[2] = {num_entity_dofs[2]},
  .num_entity_dofs[3] = {num_entity_dofs[3]},
  .tabulate_entity_dofs = tabulate_entity_dofs_{factory_name},
  .num_entity_closure_dofs[0] = {num_entity_closure_dofs[0]},
  .num_entity_closure_dofs[1] = {num_entity_closure_dofs[1]},
  .num_entity_closure_dofs[2] = {num_entity_closure_dofs[2]},
  .num_entity_closure_dofs[3] = {num_entity_closure_dofs[3]},
  .tabulate_entity_closure_dofs = tabulate_entity_closure_dofs_{factory_name},
  .num_sub_dofmaps = {num_sub_dofmaps},
  .sub_dofmaps = {sub_dofmaps}
}};

// End of code for dofmap {factory_name}
"""
