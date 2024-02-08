# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2020.
"""Code generation strings for a form."""

declaration = """
extern ufcx_form {factory_name};

// Helper used to create form using name which was given to the
// form in the UFL file.
// This helper is called in user c++ code.
//
extern ufcx_form* {name_from_uflfile};

// Helper used to create function space using function name
// i.e. name of the Python variable.
//
ufcx_function_space* functionspace_{name_from_uflfile}(const char* function_name);
"""

factory = """
// Code for form {factory_name}

{original_coefficient_position_init}
{dofmaps_init}
{finite_elements_init}
{form_integral_offsets_init}
{form_integrals_init}
{form_integral_ids_init}

{coefficient_names_init}
{constant_names_init}

ufcx_form {factory_name} =
{{

  .signature = {signature},
  .rank = {rank},
  .num_coefficients = {num_coefficients},
  .num_constants = {num_constants},
  .original_coefficient_position = {original_coefficient_position},

  .coefficient_name_map = {coefficient_names},
  .constant_name_map = {constant_names},

  .finite_elements = {finite_elements},
  .dofmaps = {dofmaps},

  .form_integrals = {form_integrals},
  .form_integral_ids = {form_integral_ids},
  .form_integral_offsets = form_integral_offsets_{factory_name}
}};

// Alias name
ufcx_form* {name_from_uflfile} = &{factory_name};

{value_shape_init}
ufcx_function_space* functionspace_{name_from_uflfile}(const char* function_name)
{{
{functionspace}
}}

// End of code for form {factory_name}
"""
