# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2018

declaration = """
extern ufc_coordinate_mapping {factory_name};

// Helper used to create coordinate map using name given to the
// UFL file.
extern ufc_coordinate_mapping* coordinate_mapping_{prefix};
"""

factory = """
// Code for coordinate mapping {factory_name}

int permute_dofs_{factory_name}(int* dof_list, const uint32_t cell_permutation)
{{
  {permute_dofs}
}}

int unpermute_dofs_{factory_name}(int* dof_list, const uint32_t cell_permutation)
{{
  {unpermute_dofs}
}}

ufc_coordinate_mapping {factory_name} =
{{
  .signature = {signature},
  .element_family = {family},
  .element_degree = {degree},
  .geometric_dimension = {geometric_dimension},
  .topological_dimension = {topological_dimension},
  .is_affine = {is_affine},
  .cell_shape = {cell_shape},
  .scalar_dofmap = {scalar_dofmap}
}};

ufc_coordinate_mapping* coordinate_mapping_{prefix} = &{factory_name};

// End of code for coordinate mapping {factory_name}
"""
