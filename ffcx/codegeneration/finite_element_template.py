# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2018.

declaration = """
extern ufc_finite_element {factory_name};
"""

factory = """
// Code for element {factory_name}

{value_shape_init}
{reference_value_shape_init}
{sub_elements_init}

int apply_dof_transformation_{factory_name}(
     double* restrict data, uint32_t cell_permutation, int dim)
{{
  {apply_dof_transformation}
}}

int apply_dof_transformation_to_scalar_{factory_name}(
     ufc_scalar_t* restrict data, uint32_t cell_permutation, int dim)
{{
  {apply_dof_transformation_to_scalar}
}}

int apply_inverse_transpose_dof_transformation_{factory_name}(
     double* restrict data, uint32_t cell_permutation, int dim)
{{
  {apply_inverse_transpose_dof_transformation}
}}

int apply_inverse_transpose_dof_transformation_to_scalar_{factory_name}(
     ufc_scalar_t* restrict data, uint32_t cell_permutation, int dim)
{{
  {apply_inverse_transpose_dof_transformation_to_scalar}
}}

ufc_finite_element {factory_name} =
{{
  .signature = {signature},
  .cell_shape = {cell_shape},
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
  .family = {family},
  .block_size = {block_size},

  .needs_transformation_data = {needs_transformation_data},
  .interpolation_is_identity = {interpolation_is_identity},

  .apply_dof_transformation = apply_dof_transformation_{factory_name},
  .apply_dof_transformation_to_scalar = apply_dof_transformation_to_scalar_{factory_name},
  .apply_inverse_transpose_dof_transformation
      = apply_inverse_transpose_dof_transformation_{factory_name},
  .apply_inverse_transpose_dof_transformation_to_scalar
      = apply_inverse_transpose_dof_transformation_to_scalar_{factory_name},
  .num_sub_elements = {num_sub_elements},
  .sub_elements = {sub_elements}
}};

// End of code for element {factory_name}
"""
