# Copyright (C) 2025 Chris Richardson
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Template for file output."""

factory = """
# Code for form {factory_name}

{original_coefficient_position_init}
{finite_element_hashes_init}
{form_integral_offsets_init}
{form_integrals_init}
{form_integral_ids_init}

{coefficient_names_init}
{constant_names_init}
{constant_ranks_init}
{constant_shapes_init}

class {factory_name}(object):

  signature = {signature}
  rank = {rank}

  num_coefficients = {num_coefficients}
  original_coefficient_positions = {original_coefficient_positions}
  coefficient_name_map = {coefficient_names}

  num_constants = {num_constants}
  constant_ranks = {constant_ranks}
  constant_shapes = {constant_shapes}
  constant_name_map = {constant_names}

  finite_element_hashes = {finite_element_hashes}

  form_integrals = {form_integrals}
  form_integral_ids = {form_integral_ids}
  form_integral_offsets = form_integral_offsets_{factory_name}


# Alias name
{name_from_uflfile} = {factory_name}

# End of code for form {factory_name}
"""
