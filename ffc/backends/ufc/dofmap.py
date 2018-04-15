# -*- coding: utf-8 -*-
# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2018.

declaration = """
ufc_dofmap* create_{factory_name}();
"""

factory = """
// Code for dofmap {factory_name}

int64_t num_entity_dofs_{factory_name}(int64_t d)
{{
{num_entity_dofs}
}}

int64_t num_entity_closure_dofs_{factory_name}(int64_t d)
{{
{num_entity_closure_dofs}
}}

void tabulate_dofs_{factory_name}(int64_t * dofs, const int64_t* num_global_entities,
                                  const int64_t** entity_indices)
{{
{tabulate_dofs}
}}

void tabulate_facet_dofs_{factory_name}(int64_t * dofs, int64_t facet)
{{
{tabulate_facet_dofs}
}}

void tabulate_entity_dofs_{factory_name}(int64_t * dofs, int64_t d, int64_t i)
{{
{tabulate_entity_dofs}
}}

void tabulate_entity_closure_dofs_{factory_name}(int64_t * dofs, int64_t d, int64_t i)
{{
{tabulate_entity_closure_dofs}
}}

ufc_dofmap* create_sub_dofmap_{factory_name}(int64_t i)
{{
{create_sub_dofmap}
}}

ufc_dofmap* create_{factory_name}()
{{
  ufc_dofmap* dofmap = malloc(sizeof(*dofmap));
  dofmap->signature = {signature};
  dofmap->num_global_support_dofs = {num_global_support_dofs};
  dofmap->num_element_support_dofs = {num_element_support_dofs};
  dofmap->num_element_dofs = {num_element_dofs};
  dofmap->num_facet_dofs = {num_facet_dofs};
  dofmap->num_entity_dofs = num_entity_dofs_{factory_name};
  dofmap->num_entity_closure_dofs = num_entity_closure_dofs_{factory_name};
  dofmap->tabulate_dofs = tabulate_dofs_{factory_name};
  dofmap->tabulate_facet_dofs = tabulate_facet_dofs_{factory_name};
  dofmap->tabulate_entity_dofs = tabulate_entity_dofs_{factory_name};
  dofmap->tabulate_entity_closure_dofs = tabulate_entity_closure_dofs_{factory_name};
  dofmap->num_sub_dofmaps = {num_sub_dofmaps};
  dofmap->create_sub_dofmap = create_sub_dofmap_{factory_name};
  dofmap->create = create_{factory_name};

  return dofmap;
}};

// End of code for dofmap {factory_name}
"""
