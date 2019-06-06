# -*- coding: utf-8 -*-
# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2018

declaration = """
ufc_integral* create_{factory_name}(void);
"""

custom_declaration = """
ufc_custom_integral* create_{factory_name}(void);
"""

tabulate_implementation = {
    "cell":
    """
void tabulate_tensor_{factory_name}(ufc_scalar_t* restrict A, const ufc_scalar_t* w,
                                    const double* restrict coordinate_dofs,
                                    const int* unused_local_index,
                                    const int* cell_orientation)
{{
{tabulate_tensor}
}}
""",
    "exterior_facet":
    """
void tabulate_tensor_{factory_name}(ufc_scalar_t* restrict A, const ufc_scalar_t* w,
                                    const double* restrict coordinate_dofs,
                                    const int* facet,
                                    const int* cell_orientation)
{{
{tabulate_tensor}
}}
""",
    "interior_facet":
    """
void tabulate_tensor_{factory_name}(ufc_scalar_t* restrict A, const ufc_scalar_t* w,
                                    const double* restrict coordinate_dofs,
                                    const int* facet,
                                    const int* cell_orientation)
{{
{tabulate_tensor}
}}
""",
    "vertex":
    """
void tabulate_tensor_{factory_name}(ufc_scalar_t* restrict A, const ufc_scalar_t* w,
                                    const double* restrict coordinate_dofs,
                                    const int* vertex,
                                    const int* cell_orientation)
{{
{tabulate_tensor}
}}
""",
    "custom":
    """
void tabulate_tensor_{factory_name}(ufc_scalar_t* restrict A, const ufc_scalar_t* w,
                          const double* restrict coordinate_dofs,
                          int num_quadrature_points,
                          const double* restrict quadrature_points,
                          const double* restrict quadrature_weights,
                          const double* restrict facet_normals,
                          int cell_orientation)
{{
{tabulate_tensor}
}}
"""
}

factory = """
// Code for integral {factory_name}

{tabulate_tensor}

ufc_integral* create_{factory_name}(void)
{{
  static const bool enabled{enabled_coefficients}

  ufc_integral* integral = malloc(sizeof(*integral));
  integral->enabled_coefficients = enabled;
  integral->tabulate_tensor = tabulate_tensor_{factory_name};
  return integral;
}};

// End of code for integral {factory_name}
"""

custom_factory = """
// Code for custom integral {factory_name}

{tabulate_tensor}

ufc_custom_integral* create_{factory_name}(void)
{{
  static const bool enabled{enabled_coefficients}

  ufc_custom_integral* integral = malloc(sizeof(*integral));
  integral->enabled_coefficients = enabled;
  integral->tabulate_tensor = tabulate_tensor_{factory_name};
  return integral;
}};

// End of code for custom integral {factory_name}
"""
