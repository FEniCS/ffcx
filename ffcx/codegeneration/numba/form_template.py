# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2020.


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

# Name: {name_from_uflfile}
# End of code for form {factory_name}
"""
