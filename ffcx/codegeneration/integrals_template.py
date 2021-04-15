# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2018

declaration = """
extern ufc_integral {factory_name};
"""

factory = """
// Code for integral {factory_name}

void tabulate_tensor_{factory_name}(ufc_scalar_t* restrict A,
                                    const ufc_scalar_t* restrict w,
                                    const ufc_scalar_t* restrict c,
                                    const double* restrict coordinate_dofs,
                                    const int* restrict entity_local_index,
                                    const uint8_t* restrict quadrature_permutation,
                                    const uint32_t cell_permutation)
{{
{tabulate_tensor}
}}

{enabled_coefficients_init}

ufc_integral {factory_name} =
{{
  .enabled_coefficients = {enabled_coefficients},
  .tabulate_tensor = tabulate_tensor_{factory_name},
  .needs_transformation_data = {needs_transformation_data}
}};

// End of code for integral {factory_name}
"""
