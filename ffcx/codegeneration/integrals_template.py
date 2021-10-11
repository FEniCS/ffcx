# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2018

declaration = """
extern ufc_integral {factory_name};
"""

factory = """
// Code for integral {factory_name}

void tabulate_tensor_{factory_name}({scalar_type}* restrict A,
                                    const {scalar_type}* restrict w,
                                    const {scalar_type}* restrict c,
                                    const double* restrict coordinate_dofs,
                                    const int* restrict entity_local_index,
                                    const uint8_t* restrict quadrature_permutation)
{{
{tabulate_tensor}
}}

{enabled_coefficients_init}

ufc_integral {factory_name} =
{{
  .enabled_coefficients = {enabled_coefficients},
  .tabulate_tensor = tabulate_tensor_{factory_name},
  .needs_facet_permutations = {needs_facet_permutations},
}};

// End of code for integral {factory_name}
"""
