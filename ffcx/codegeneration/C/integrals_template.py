# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2018
"""Code generation strings for an integral."""

declaration = """
extern ufcx_integral {factory_name};
"""

factory = """
// Code for integral {factory_name}

void tabulate_tensor_{factory_name}({scalar_type}* restrict A,
                                    const {scalar_type}* restrict w,
                                    const {scalar_type}* restrict c,
                                    const {geom_type}* restrict coordinate_dofs,
                                    const int* restrict entity_local_index,
                                    const uint8_t* restrict quadrature_permutation)
{{
{tabulate_tensor}
}}

{enabled_coefficients_init}

ufcx_integral {factory_name} =
{{
  .enabled_coefficients = {enabled_coefficients},
  .tabulate_tensor_float32 = {tabulate_tensor_float32},
  .tabulate_tensor_float64 = {tabulate_tensor_float64},
  .tabulate_tensor_complex64 = {tabulate_tensor_complex64},
  .tabulate_tensor_complex128 = {tabulate_tensor_complex128},
  .needs_facet_permutations = {needs_facet_permutations},
  .coordinate_element = {coordinate_element},
}};

// End of code for integral {factory_name}
"""
