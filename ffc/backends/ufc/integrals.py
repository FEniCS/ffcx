# -*- coding: utf-8 -*-
# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2018

ufc_integral_declaration = """
ufc_{type}_integral* create_{factory_name}();
"""

tabulate_tensor_implementation = {
    "cell":
    """
void tabulate_tensor_{factory_name}(double *A, const double *const *w,
                                    const double *coordinate_dofs,
                                    int cell_orientation)
{{
{tabulate_tensor}
}}
""",
    "exterior_facet":
    """
void tabulate_tensor_{factory_name}(double *A, const double *const *w,
                                     const double *coordinate_dofs,
                                     int64_t facet,
                                     int cell_orientation)
{{
{tabulate_tensor}
}}
""",
    "interior_facet":
    """
void tabulate_tensor_{factory_name}(double *A, const double *const *w,
                          const double *coordinate_dofs_0,
                          const double *coordinate_dofs_1, int64_t facet_0,
                          int64_t facet_1, int cell_orientation_0,
                          int cell_orientation_1)
{{
{tabulate_tensor}
}}
""",
    "vertex":
    """
void tabulate_tensor_{factory_name}(double *A, const double *const *w,
                          const double *coordinate_dofs, int64_t vertex,
                          int cell_orientation)
{{
{tabulate_tensor}
}}
""",
    "custom":
    """
void tabulate_tensor_{factory_name}(double *A, const double *const *w,
                          const double *coordinate_dofs,
                          int64_t num_quadrature_points,
                          const double *quadrature_points,
                          const double *quadrature_weights,
                          const double *facet_normals,
                          int cell_orientation)
{{
{tabulate_tensor}
}}
"""
}

ufc_integral_factory = """
// Code for {type}_integral {factory_name}

{tabulate_tensor}

ufc_{type}_integral* create_{factory_name}()
{{
  static const bool enabled{enabled_coefficients}

  ufc_{type}_integral* integral = malloc(sizeof(*integral));
  integral->enabled_coefficients = enabled;
  integral->tabulate_tensor = tabulate_tensor_{factory_name};
  return integral;
}};

// End of code for {type}_integral {factory_name}
"""
