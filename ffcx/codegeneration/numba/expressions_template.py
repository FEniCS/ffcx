# Copyright (C) 2019 Michal Habera
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Template for expression output."""

declaration = """
"""

factory = """
# Code for expression {factory_name}

def tabulate_tensor_{factory_name}(_A, _w, _c, _coordinate_dofs,
                                   _entity_local_index,
                                   _quadrature_permutation, custom_data):
{tabulate_expression}



class {factory_name}:

    tabulate_tensor = tabulate_tensor_{factory_name}
    num_coefficients = {num_coefficients}
    num_constants = {num_constants}
    original_coefficient_positions = {original_coefficient_positions}
    coefficient_names = {coefficient_names}
    constant_names = {constant_names}
    num_points = {num_points}
    topological_dimension = {topological_dimension}
    points = {points}
    value_shape = {value_shape}
    num_components = {num_components}
    rank = {rank}
    # function_spaces = {function_spaces}

{name_from_uflfile} = {factory_name}

# Name: {name_from_uflfile}
# End of code for expression {factory_name}
"""
