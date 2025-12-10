# Copyright (C) 2025 Chris Richardson and Paul T. Kühner
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Templates for C++ integral output."""

declaration = """
class {factory_name}
{{
public:

// Constructor
{factory_name}();

// Kernel
template <typename T, typename U>
void tabulate_tensor(T* A,
                     const T* RESTRICT w,
                     const T* RESTRICT c,
                     const U* RESTRICT coordinate_dofs,
                     const std::int32_t* RESTRICT entity_local_index,
                     const std::uint8_t* RESTRICT quadrature_permutation);

// Data
std::vector<bool> enabled_coefficients; // TODO: std::vector<char>?
bool needs_facet_permutations;

}};
"""

factory = """
// Code for integral {factory_name}

template <typename T, typename U>
void {factory_name}::tabulate_tensor(T* RESTRICT A,
                     const T* RESTRICT w,
                     const T* RESTRICT c,
                     const U* RESTRICT coordinate_dofs,
                     const std::int32_t* RESTRICT entity_local_index,
                     const std::uint8_t* RESTRICT quadrature_permutation)
{{
{tabulate_tensor}
}}

{factory_name}::{factory_name}()
{{
  enabled_coefficients = {enabled_coefficients};
  needs_facet_permutations = {needs_facet_permutations};
}}

// End of code for integral {factory_name}
"""
