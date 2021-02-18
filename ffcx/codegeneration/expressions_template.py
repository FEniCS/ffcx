# Copyright (C) 2019 Michal Habera
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

declaration = """
ufc_expression* create_{factory_name}(void);

// Helper used to create expression
// This helper is called in user c++ code.
//
ufc_expression* create_expression(void);

"""

factory = """
// Code for expression {factory_name}

void tabulate_tensor_{factory_name}(ufc_scalar_t* restrict A,
                                    const ufc_scalar_t* restrict w,
                                    const ufc_scalar_t* restrict c,
                                    const double* restrict coordinate_dofs)
{{
{tabulate_expression}
}}

ufc_expression* create_{factory_name}(void)
{{

  ufc_expression* expression = (ufc_expression*)malloc(sizeof(*expression));

{original_coefficient_positions}
{coefficient_names}
{constant_names}
{points}
{value_shape}

  expression->tabulate_expression = tabulate_tensor_{factory_name};
  expression->num_coefficients = {num_coefficients};
  expression->num_points = {num_points};
  expression->topological_dimension = {topological_dimension};
  expression->points = *points;
  expression->value_shape = value_shape;
  expression->num_components = {num_components};

  return expression;
}}

ufc_expression* create_expression(void)
{{
  return create_{factory_name}();
}}

// End of code for expression {factory_name}
"""
