# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2018.

factory = """
# Code for dofmap {factory_name}

class {factory_name}(object):
    signature = {signature}
    num_global_support_dofs = {num_global_support_dofs}
    num_element_support_dofs = {num_element_support_dofs}
    block_size = {block_size}
    entity_dofs = {entity_dofs}
    entity_closure_dofs = {entity_closure_dofs}
    num_sub_dofmaps = {num_sub_dofmaps}
    sub_dofmaps = {sub_dofmaps}

# End of code for dofmap {factory_name}
"""
