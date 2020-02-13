# Copyright (C) 2019 Michal Habera
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

declaration = """
ufc_expression* create_expression_{factory_name}(void);
"""

factory = """
// Code for expression {factory_name}

void tabulate_expression_{factory_name}(ufc_scalar_t* restrict A, const ufc_scalar_t* w,
                                        const ufc_scalar_t* c,
                                        const double* restrict coordinate_dofs)
{{
{tabulate_expression}
}}

ufc_expression* create_{factory_name}(void)
{{

  {original_coefficient_positions}
  {points}
  {value_shape}

  ufc_expression* expression = malloc(sizeof(*expression));
  expression->tabulate_expression = tabulate_expression_{factory_name};
  expression->original_coefficient_positions = original_coefficient_positions;
  expression->num_coefficients = {num_coefficients};
  expression->num_points = {num_points};
  expression->topological_dimension = {topological_dimension};
  expression->points = *points;
  expression->value_shape = value_shape;
  expression->num_components = {num_components};

  return expression;
}}

// End of code for expression {factory_name}
"""
