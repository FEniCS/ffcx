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

int transform_reference_basis_derivatives_{factory_name}(
    double * restrict values, int order, int num_points,
    const double * restrict reference_values,
    const double * restrict X, const double * restrict J,
    const double * restrict detJ, const double * restrict K,
    const uint32_t cell_permutation)
{{
  {transform_reference_basis_derivatives}
}}

int transform_values_{factory_name}(
     ufc_scalar_t* restrict reference_values,
     const ufc_scalar_t* restrict physical_values,
     const double* restrict coordinate_dofs,
     const ufc_coordinate_mapping* restrict cm)
{{
  {transform_values}
}}

int interpolate_into_cell_{factory_name}(ufc_scalar_t* restrict coefficients, const ufc_scalar_t* restrict evaluations,
      const uint32_t cell_permutation)
{{
  {interpolate_into_cell}
}}

int tabulate_reference_dof_coordinates_{factory_name}(double* restrict reference_dof_coordinates)
{{
  {tabulate_reference_dof_coordinates}
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

  element->needs_permutation_data = {needs_permutation_data};

  element->num_interpolation_points = {num_interpolation_points};
  {interpolation_points};
  element->interpolate_into_cell = interpolate_into_cell_{factory_name};

  element->transform_reference_basis_derivatives = transform_reference_basis_derivatives_{factory_name};
  element->transform_values = transform_values_{factory_name};
  element->tabulate_reference_dof_coordinates = tabulate_reference_dof_coordinates_{factory_name};
  element->num_sub_elements = {num_sub_elements};
  element->create_sub_element = create_sub_element_{factory_name};
  element->create = create_{factory_name};

  return element;
}}

// End of code for element {factory_name}
"""
