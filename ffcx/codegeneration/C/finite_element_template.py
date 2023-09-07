# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2018.

declaration = """
extern ufcx_finite_element {factory_name};
"""

factory = """
// Code for element {factory_name}

{value_shape_init}
{reference_value_shape_init}
{sub_elements_init}
{custom_element_init}

ufcx_finite_element {factory_name} =
{{
  .signature = {signature},
  .cell_shape = {cell_shape},
  .element_type = {element_type},
  .topological_dimension = {topological_dimension},
  .geometric_dimension = {geometric_dimension},
  .space_dimension = {space_dimension},
  .value_rank = {value_rank},
  .value_shape = {value_shape},
  .value_size = {value_size},
  .reference_value_rank = {reference_value_rank},
  .reference_value_shape = {reference_value_shape},
  .reference_value_size = {reference_value_size},
  .degree = {degree},
  .block_size = {block_size},
  .family = {family},
  .basix_family = {basix_family},
  .basix_cell = {basix_cell},
  .discontinuous = {discontinuous},
  .lagrange_variant = {lagrange_variant},
  .dpc_variant = {dpc_variant},
  .num_sub_elements = {num_sub_elements},
  .sub_elements = {sub_elements},
  .custom_element = {custom_element}
}};

// End of code for element {factory_name}
"""
