# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2018.

declaration = """
ufc_finite_element* create_{factory_name}(void);
"""

factory = """
// Code for element {factory_name}

int value_dimension_{factory_name}(int i)
{{
  {value_dimension}
}}

int reference_value_dimension_{factory_name}(int i)
{{
  {reference_value_dimension}
}}

{sub_element_declaration}
ufc_finite_element* create_sub_element_{factory_name}(int i)
{{
  {create_sub_element}
}}

ufc_finite_element* create_{factory_name}(void)
{{
  ufc_finite_element* element = (ufc_finite_element*)malloc(sizeof(*element));

  element->signature = {signature};
  element->cell_shape = {cell_shape};
  element->topological_dimension = {topological_dimension};
  element->geometric_dimension = {geometric_dimension};
  element->space_dimension = {space_dimension};
  element->value_rank = {value_rank};
  element->value_dimension = value_dimension_{factory_name};
  element->value_size = {value_size};
  element->reference_value_rank = {reference_value_rank};
  element->reference_value_dimension = reference_value_dimension_{factory_name};
  element->reference_value_size = {reference_value_size};
  element->degree = {degree};
  element->family = {family};
  element->block_size = {block_size};

  element->needs_transformation_data = {needs_transformation_data};
  element->interpolation_is_identity = {interpolation_is_identity};

  element->num_sub_elements = {num_sub_elements};
  element->create_sub_element = create_sub_element_{factory_name};
  element->create = create_{factory_name};

  return element;
}}

// End of code for element {factory_name}
"""
