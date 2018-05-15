# -*- coding: utf-8 -*-
# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2018.

declaration = """
ufc_form* create_{factory_name}(void);
"""

factory = """
// Code for form {factory_name}

int original_coefficient_position_{factory_name}(int i)
{{
{original_coefficient_position}
}}

ufc_finite_element* create_coordinate_finite_element_{factory_name}(void)
{{
{create_coordinate_finite_element}
}}

ufc_dofmap* create_coordinate_dofmap_{factory_name}(void)
{{
{create_coordinate_dofmap}
}}

ufc_coordinate_mapping* create_coordinate_mapping_{factory_name}(void)
{{
{create_coordinate_mapping}
}}

ufc_finite_element* create_finite_element_{factory_name}(int i)
{{
{create_finite_element}
}}

ufc_dofmap* create_dofmap_{factory_name}(int i)
{{
{create_dofmap}
}}

ufc_cell_integral* create_cell_integral_{factory_name}(int subdomain_id)
{{
  {create_cell_integral}
}}

ufc_exterior_facet_integral* create_exterior_facet_integral_{factory_name}(int subdomain_id)
{{
  {create_exterior_facet_integral}
}}

ufc_interior_facet_integral* create_interior_facet_integral_{factory_name}(int subdomain_id)
{{
{create_interior_facet_integral}
}}

ufc_vertex_integral* create_vertex_integral_{factory_name}(int subdomain_id)
{{
{create_vertex_integral}
}}

ufc_custom_integral* create_custom_integral_{factory_name}(int subdomain_id)
{{
{create_custom_integral}
}}

ufc_cell_integral* create_default_cell_integral_{factory_name}(void)
{{
{create_default_cell_integral}
}}

ufc_exterior_facet_integral* create_default_exterior_facet_integral_{factory_name}(void)
{{
{create_default_exterior_facet_integral}
}}

ufc_interior_facet_integral* create_default_interior_facet_integral_{factory_name}(void)
{{
{create_default_interior_facet_integral}
}}

ufc_vertex_integral* create_default_vertex_integral_{factory_name}(void)
{{
{create_default_vertex_integral}
}}

ufc_custom_integral* create_default_custom_integral_{factory_name}(void)
{{
{create_default_custom_integral}
}}

ufc_form* create_{factory_name}(void)
{{
  ufc_form* form = malloc(sizeof(*form));

  form->signature = {signature};
  form->rank = {rank};
  form->num_coefficients = {num_coefficients};
  form->original_coefficient_position = original_coefficient_position_{factory_name};
  form->create_coordinate_finite_element = create_coordinate_finite_element_{factory_name};
  form->create_coordinate_dofmap = create_coordinate_dofmap_{factory_name};
  form->create_coordinate_mapping = create_coordinate_mapping_{factory_name};
  form->create_finite_element = create_finite_element_{factory_name};
  form->create_dofmap = create_dofmap_{factory_name};

  form->max_cell_subdomain_id = {max_cell_subdomain_id};
  form->max_exterior_facet_subdomain_id = {max_exterior_facet_subdomain_id};
  form->max_interior_facet_subdomain_id = {max_interior_facet_subdomain_id};
  form->max_vertex_subdomain_id = {max_vertex_subdomain_id};
  form->max_custom_subdomain_id = {max_custom_subdomain_id};
  form->has_cell_integrals = {has_cell_integrals};
  form->has_exterior_facet_integrals = {has_exterior_facet_integrals};
  form->has_interior_facet_integrals = {has_interior_facet_integrals};
  form->has_vertex_integrals = {has_vertex_integrals};
  form->has_custom_integrals = {has_custom_integrals};

  form->create_cell_integral = create_cell_integral_{factory_name};
  form->create_exterior_facet_integral = create_exterior_facet_integral_{factory_name};
  form->create_interior_facet_integral = create_interior_facet_integral_{factory_name};
  form->create_vertex_integral = create_vertex_integral_{factory_name};
  form->create_custom_integral = create_custom_integral_{factory_name};

  form->create_default_cell_integral = create_default_cell_integral_{factory_name};
  form->create_default_exterior_facet_integral = create_default_exterior_facet_integral_{factory_name};
  form->create_default_interior_facet_integral = create_default_interior_facet_integral_{factory_name};
  form->create_default_vertex_integral = create_default_vertex_integral_{factory_name};
  form->create_default_custom_integral = create_default_custom_integral_{factory_name};

  return form;
}};

// End of code for form {factory_name}
"""
