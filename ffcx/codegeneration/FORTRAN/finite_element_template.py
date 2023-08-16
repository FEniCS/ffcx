# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2018.

declaration = """
extern ufcx_finite_element {factory_name};
"""

factory = """
// Code for element {factory_name}

{sub_elements_init}
{custom_element_init}

type(ufcx_finite_element) :: {factory_name}

{factory_name}%signature = {signature}
{factory_name}%cell_shape = {cell_shape}
{factory_name}%element_type = {element_type}
{factory_name}%topological_dimension = {topological_dimension}
{factory_name}%geometric_dimension = {geometric_dimension}
{factory_name}%space_dimension = {space_dimension}
{factory_name}%value_rank = {value_rank}
{factory_name}%value_shape = {value_shape}
{factory_name}%value_size = {value_size}
{factory_name}%reference_value_rank = {reference_value_rank}
{factory_name}%reference_value_shape = {reference_value_shape}
{factory_name}%reference_value_size = {reference_value_size}
{factory_name}%degree = {degree}
{factory_name}%block_size = {block_size}
{factory_name}%family = {family}
{factory_name}%basix_family = {basix_family}
{factory_name}%basix_cell = {basix_cell}
{factory_name}%discontinuous = {discontinuous}
{factory_name}%lagrange_variant = {lagrange_variant}
{factory_name}%dpc_variant = {dpc_variant}
{factory_name}%num_sub_elements = {num_sub_elements}
{factory_name}%sub_elements = {sub_elements}
{factory_name}%custom_element = {custom_element}


// End of code for element {factory_name}
"""
