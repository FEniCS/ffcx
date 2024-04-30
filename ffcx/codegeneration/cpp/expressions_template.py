# Copyright (C) 2019 Michal Habera
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

declaration = """
extern ufcx_expression {factory_name};

// Helper used to create expression using name which was given to the
// expression in the UFL file.
// This helper is called in user c++ code.
//
extern ufcx_expression* {name_from_uflfile};
"""

factory = """
// Code for expression {factory_name}

void tabulate_tensor_{factory_name}({scalar_type}* restrict A,
                                    const {scalar_type}* restrict w,
                                    const {scalar_type}* restrict c,
                                    const {geom_type}* restrict coordinate_dofs,
                                    const int* restrict entity_local_index,
                                    const uint8_t* restrict quadrature_permutation)
{{
{tabulate_expression}
}}

{points_init}
{value_shape_init}
{original_coefficient_positions_init}
{function_spaces_alloc}
{function_spaces_init}
{coefficient_names_init}
{constant_names_init}


ufcx_expression {factory_name} =
{{
  .tabulate_tensor_{np_scalar_type} = tabulate_tensor_{factory_name},
  .num_coefficients = {num_coefficients},
  .num_constants = {num_constants},
  .original_coefficient_positions = {original_coefficient_positions},
  .coefficient_names = {coefficient_names},
  .constant_names = {constant_names},
  .num_points = {num_points},
  .topological_dimension = {topological_dimension},
  .points = {points},
  .value_shape = {value_shape},
  .num_components = {num_components},
  .rank = {rank},
  .function_spaces = {function_spaces}
}};

// Alias name
ufcx_expression* {name_from_uflfile} = &{factory_name};

// End of code for expression {factory_name}
"""
