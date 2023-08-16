# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2020.

declaration = """
type(ufcx_form) :: {factory_name}
"""

factory = """
! Code for form {factory_name}

{factory_name} = ufcx_form({signature}, {rank}, {num_coefficients},
      {num_constants},
      {original_coefficient_position},
      coefficient_name_{factory_name},
      constant_name_{factory_name},
      {finite_elements},
      {dofmaps},
      integral_ids_{factory_name},
      num_integrals_{factory_name},
      integrals_{factory_name})


! End of code for form {factory_name}
"""
