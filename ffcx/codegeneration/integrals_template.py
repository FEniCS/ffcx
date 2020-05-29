# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2018

declaration = """
int init_{factory_name}(ufc_integral* integral);
void destroy_{factory_name}(ufc_integral* integral);
ufc_integral* create_{factory_name}(void);
"""

custom_declaration = """
int init_{factory_name}(ufc_custom_integral* integral);
void destroy_{factory_name}(ufc_custom_integral* integral);
ufc_custom_integral* create_{factory_name}(void);
"""

tabulate_implementation = {
    "cell":
    """
void tabulate_tensor_{factory_name}(ufc_scalar_t* restrict A,
                                    const ufc_scalar_t* restrict w,
                                    const ufc_scalar_t* restrict c,
                                    const double* restrict coordinate_dofs,
                                    const int* restrict unused_local_index,
                                    const uint8_t* restrict quadrature_permutation,
                                    const uint32_t cell_permutation)
{{
{tabulate_tensor}
}}
""",
    "exterior_facet":
    """
void tabulate_tensor_{factory_name}(ufc_scalar_t* restrict A,
                                    const ufc_scalar_t* restrict w,
                                    const ufc_scalar_t* restrict c,
                                    const double* restrict coordinate_dofs,
                                    const int* restrict facet,
                                    const uint8_t* restrict quadrature_permutation,
                                    const uint32_t cell_permutation)
{{
{tabulate_tensor}
}}
""",
    "interior_facet":
    """
void tabulate_tensor_{factory_name}(ufc_scalar_t* restrict A,
                                    const ufc_scalar_t* restrict w,
                                    const ufc_scalar_t* restrict c,
                                    const double* restrict coordinate_dofs,
                                    const int* restrict facet,
                                    const uint8_t* restrict quadrature_permutation,
                                    const uint32_t cell_permutation)
{{
{tabulate_tensor}
}}
""",
    "vertex":
    """
void tabulate_tensor_{factory_name}(ufc_scalar_t* restrict A,
                                    const ufc_scalar_t* restrict w,
                                    const ufc_scalar_t* restrict c,
                                    const double* restrict coordinate_dofs,
                                    const int* restrict vertex,
                                    const uint8_t* restrict quadrature_permutation,
                                    const uint32_t cell_permutation)
{{
{tabulate_tensor}
}}
""",
    "custom":
    """
void tabulate_tensor_{factory_name}(ufc_scalar_t* restrict A,
                                    const ufc_scalar_t* restrict w,
                                    const ufc_scalar_t* restrict c,
                                    const double* restrict coordinate_dofs,
                                    int num_quadrature_points,
                                    const double* restrict quadrature_points,
                                    const double* restrict quadrature_weights,
                                    const double* restrict facet_normals)
{{
{tabulate_tensor}
}}
"""
}

factory = """
// Code for integral {factory_name}

{tabulate_tensor}

void destroy_{factory_name}(ufc_integral* integral);
ufc_integral* create_{factory_name}(void);

int init_{factory_name}(ufc_integral* integral)
{{
  static const bool enabled{enabled_coefficients}
  integral->enabled_coefficients = enabled;
  integral->tabulate_tensor = tabulate_tensor_{factory_name};
  integral->init = init_{factory_name};
  integral->destroy = destroy_{factory_name};
  integral->create = create_{factory_name};
  return 0;
}}

void destroy_{factory_name}(ufc_integral* integral)
{{
}}

ufc_integral* create_{factory_name}(void)
{{
  ufc_integral* integral = malloc(sizeof(*integral));
  init_{factory_name}(integral);
  return integral;
}}

// End of code for integral {factory_name}
"""

custom_factory = """
// Code for custom integral {factory_name}

{tabulate_tensor}

void destroy_{factory_name}(ufc_custom_integral* integral);
ufc_custom_integral* create_{factory_name}(void);

int init_{factory_name}(ufc_custom_integral* integral)
{{
  static const bool enabled{enabled_coefficients}
  integral->enabled_coefficients = enabled;
  integral->tabulate_tensor = tabulate_tensor_{factory_name};
  integral->init = init_{factory_name};
  integral->destroy = destroy_{factory_name};
  integral->create = create_{factory_name};
  return 0;
}}

void destroy_{factory_name}(ufc_custom_integral* integral)
{{
}}

ufc_custom_integral* create_{factory_name}(void)
{{
  ufc_custom_integral* integral = malloc(sizeof(*integral));
  init_{factory_name}(integral);
  return integral;
}}

// End of code for custom integral {factory_name}
"""
