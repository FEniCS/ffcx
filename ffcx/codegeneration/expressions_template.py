# Copyright (C) 2019 Michal Habera
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

declaration = """
extern ufcx_expression {factory_name};
"""

factory = """
// Code for expression {factory_name}

void tabulate_expression_{factory_name}({scalar_type}* restrict A,
                                        const {scalar_type}* restrict w,
                                        const {scalar_type}* restrict c,
                                        const double* restrict coordinate_dofs)
{{
{tabulate_expression}
}}

{points_init}
{value_shape_init}
{original_coefficient_positions_init}

ufcx_expression {factory_name} =
{{
  .tabulate_expression = tabulate_expression_{factory_name},
  .num_coefficients = {num_coefficients},
  .num_points = {num_points},
  .topological_dimension = {topological_dimension},
  .needs_facet_permutations = {needs_facet_permutations},
  .points = {points},
  .value_shape = {value_shape},
  .num_components = {num_components},
  .original_coefficient_positions = {original_coefficient_positions}
}};

// End of code for expression {factory_name}
"""
