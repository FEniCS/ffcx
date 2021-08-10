# Copyright (C) 2019 Michal Habera
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

declaration = """
<<<<<<< HEAD
ufc_expression* create_{factory_name}(void);

// Helper used to create expression
// This helper is called in user c++ code.
//
ufc_expression* create_expression(void);

=======
extern ufc_expression {factory_name};
>>>>>>> main
"""

factory = """
// Code for expression {factory_name}

void tabulate_tensor_{factory_name}(ufc_scalar_t* restrict A,
                                    const ufc_scalar_t* restrict w,
                                    const ufc_scalar_t* restrict c,
                                    const double* restrict coordinate_dofs,
                                    const int* restrict unused_local_index,
                                    const uint8_t* restrict quadrature_permutation,
                                    const uint32_t cell_permutation)
{{
{tabulate_expression}
}}

{points_init}
{value_shape_init}
{original_coefficient_positions_init}

ufc_expression {factory_name} =
{{
  .tabulate_expression = tabulate_expression_{factory_name},
  .num_coefficients = {num_coefficients},
  .num_points = {num_points},
  .topological_dimension = {topological_dimension},
  .points = {points},
  .value_shape = {value_shape},
  .num_components = {num_components},
  .original_coefficient_positions = {original_coefficient_positions}
}};

ufc_expression* create_expression(void)
{{
  return create_{factory_name}();
}}

// End of code for expression {factory_name}
"""
