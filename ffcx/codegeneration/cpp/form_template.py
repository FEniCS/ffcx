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

{dofmaps_init}
{finite_elements_init}

{name_from_uflfile}::constant_name = {constant_name_map};
{name_from_uflfile}::coefficient_name = {coefficient_name_map};
{name_from_uflfile}::signature ={signature};
{name_from_uflfile}::rank = {rank};
{name_from_uflfile}::num_coefficients = {num_coefficients};
{name_from_uflfile}::num_constants = {num_constants};
{name_from_uflfile}::original_coefficient_position = {original_coefficient_position};
{name_from_uflfile}::coefficient_name_map = coefficient_name_{factory_name};
{name_from_uflfile}::constant_name_map = constant_name_{factory_name};
{name_from_uflfile}::finite_elements = {finite_elements};
{name_from_uflfile}::dofmaps = {dofmaps};
{name_from_uflfile}::form_integrals = {form_integrals};
{name_from_uflfile}::form_integral_ids = {form_integral_ids};
{name_from_uflfile}::form_integral_offsets = {form_integral_offsets};

ufcx_function_space* functionspace_{name_from_uflfile}(const char* function_name)
{{
{functionspace}
}}

// End of code for form {factory_name}
"""
