# -*- coding: utf-8 -*-
# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2018.

declaration = """
ufc_dofmap* create_{factory_name}(void);
"""

factory = """
// Code for dofmap {factory_name}

void tabulate_entity_dofs_{factory_name}(int* restrict dofs, int d, int i)
{{
{tabulate_entity_dofs}
}}

{sub_dofmap_declaration}
ufc_dofmap* create_sub_dofmap_{factory_name}(int i)
{{
{create_sub_dofmap}
}}

ufc_dofmap* create_{factory_name}(void)
{{
  ufc_dofmap* dofmap = malloc(sizeof(*dofmap));
  dofmap->signature = {signature};
  dofmap->num_global_support_dofs = {num_global_support_dofs};
  dofmap->num_element_support_dofs = {num_element_support_dofs};
  dofmap->num_entity_dofs[0] = {num_entity_dofs[0]};
  dofmap->num_entity_dofs[1] = {num_entity_dofs[1]};
  dofmap->num_entity_dofs[2] = {num_entity_dofs[2]};
  dofmap->num_entity_dofs[3] = {num_entity_dofs[3]};
  dofmap->tabulate_entity_dofs = tabulate_entity_dofs_{factory_name};
  dofmap->num_sub_dofmaps = {num_sub_dofmaps};
  dofmap->create_sub_dofmap = create_sub_dofmap_{factory_name};
  dofmap->create = create_{factory_name};

  return dofmap;
}};

// End of code for dofmap {factory_name}
"""
