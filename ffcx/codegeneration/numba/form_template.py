# Copyright (C) 2025 Chris Richardson
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Template for file output."""

factory = """
# Code for form {factory_name}

class {factory_name}(object):

  signature = {signature}
  rank = {rank}
  num_coefficients = {num_coefficients}
  num_constants = {num_constants}
  original_coefficient_position = {original_coefficient_position}

  coefficient_name_map = {coefficient_name_map}
  constant_name_map = {constant_name_map}

  form_integrals = {form_integrals}
  form_integral_ids = {form_integral_ids}
  form_integral_offsets = {form_integral_offsets}

{name_from_uflfile} = {factory_name}

# End of code for form {factory_name}
"""
