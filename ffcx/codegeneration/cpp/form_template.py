# Copyright (C) 2025 Chris Richardson and Paul T. Kühner
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Templates for C++ form output."""

declaration = """
extern ufcx_form {factory_name};

// Helper used to create form using name which was given to the
// form in the UFL file.
// This helper is called in user c++ code.
//
extern ufcx_form* {name_from_uflfile};

"""

factory = """
// Code for form {factory_name}

// TODO: that correct?
{original_coefficient_position_init}
{finite_element_hashes_init}
{form_integral_offsets_init}
{form_integrals_init}
{form_integral_ids_init}

{coefficient_names_init}
{constant_names_init}
{constant_ranks_init}
{constant_shapes_init}

{name_from_uflfile}::signature ={signature};
{name_from_uflfile}::rank = {rank};

{name_from_uflfile}::num_coefficients = {num_coefficients};
{name_from_uflfile}::original_coefficient_positions = {original_coefficient_positions};
{name_from_uflfile}::coefficient_name_map = {coefficient_names};

{name_from_uflfile}::num_constants = {num_constants};
{name_from_uflfile}::constant_ranks = {constant_ranks};
{name_from_uflfile}::constant_shapes = {constant_shapes};
{name_from_uflfile}::constant_name_map = {constant_names};

{name_from_uflfile}::finite_element_hashes = {finite_element_hashes},

{name_from_uflfile}::form_integrals = {form_integrals};
{name_from_uflfile}::form_integral_ids = {form_integral_ids};
{name_from_uflfile}::form_integral_offsets = form_integral_offsets_{factory_name};

// Alias name
using {name_from_uflfile} = {factory_name};

// End of code for form {factory_name}
"""
