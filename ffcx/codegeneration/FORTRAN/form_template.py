# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2020.

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

// Return a list of the coefficient names.
const char** coefficient_name_{factory_name}(void)
{{
{coefficient_name_map}
}}

// Return a list of the constant names.
const char** constant_name_{factory_name}(void)
{{
{constant_name_map}
}}

int* integral_ids_{factory_name}(ufcx_integral_type integral_type)
{{
{integral_ids}
}}

int num_integrals_{factory_name}(ufcx_integral_type integral_type)
{{
{num_integrals}
}}

ufcx_integral** integrals_{factory_name}(ufcx_integral_type integral_type)
{{
{integrals}
}}

ufcx_form {factory_name} =
{{

  .signature = {signature},
  .rank = {rank},
  .num_coefficients = {num_coefficients},
  .num_constants = {num_constants},
  .original_coefficient_position = {original_coefficient_position},

  .coefficient_name_map = coefficient_name_{factory_name},
  .constant_name_map = constant_name_{factory_name},

  .finite_elements = {finite_elements},
  .dofmaps = {dofmaps},

  .integral_ids = integral_ids_{factory_name},
  .num_integrals = num_integrals_{factory_name},

  .integrals = integrals_{factory_name}
}};

// Alias name
ufcx_form* {name_from_uflfile} = &{factory_name};

ufcx_function_space* functionspace_{name_from_uflfile}(const char* function_name)
{{
{functionspace}
}}

// End of code for form {factory_name}
"""
