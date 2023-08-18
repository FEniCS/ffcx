# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2018

declaration = """
class {factory_name} public ufcx_integral
{{
public:

template <typename T, typename U>
void tabulate_tensor(T* A,
                     const T* w,
                     const T* c,
                     const U* coordinate_dofs,
                     const int* entity_local_index,
                     const uint8_t* quadrature_permutation);

std::vector<bool> enabled_coefficients;

}};
"""

factory = """
// Code for integral {factory_name}

template <typename T, typename U>
void {factory_name}::tabulate_tensor(T* A,
                     const T* w,
                     const T* c,
                     const U* coordinate_dofs,
                     const int* entity_local_index,
                     const uint8_t* quadrature_permutation)
{{
{tabulate_tensor}
}}

{factory_name}::enabled_coefficients = {enabled_coefficients};
{factory_name}::needs_facet_permutations = {needs_facet_permutations};
{factory_name}::coordinate_element = &{coordinate_element};

// End of code for integral {factory_name}
"""
