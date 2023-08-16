# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2018.

declaration = """
type(ufcx_dofmap) :: {factory_name}
"""

factory = """
! Code for dofmap {factory_name}

{factory_name}%signature = {signature}
{factory_name}%num_global_support_dofs = {num_global_support_dofs}
{factory_name}%num_element_support_dofs = {num_element_support_dofs}
{factory_name}%block_size = {block_size}
{factory_name}%entity_dofs = {entity_dofs}
{factory_name}%entity_dof_offsets = {entity_dof_offsets}
{factory_name}%entity_closure_dofs = {entity_closure_dofs}
{factory_name}%entity_closure_dof_offsets = {entity_closure_dof_offsets}
{factory_name}%num_sub_dofmaps = {num_sub_dofmaps}
{factory_name}%sub_dofmaps = {sub_dofmaps}
! End of code for dofmap {factory_name}
"""
