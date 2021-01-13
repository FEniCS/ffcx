# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2018

declaration = """
ufc_coordinate_mapping* create_{factory_name}(void);

// Helper used to create coordinate map using name given to the
// UFL file.
// This helper is called in user c++ code.
//
ufc_coordinate_mapping* create_coordinate_map_{prefix}(void);
"""

# declaration = """
# ufc_coordinate_mapping* create_{factory_name}(void);
# """

factory = """
// Code for coordinate mapping {factory_name}

int get_dof_permutation_{factory_name}(int* dof_list, const uint32_t cell_permutation)
{{
  {get_dof_permutation}
}}

ufc_coordinate_mapping* create_{factory_name}(void)
{{
  ufc_coordinate_mapping* cmap = (ufc_coordinate_mapping*)malloc(sizeof(*cmap));
  cmap->signature = {signature};
  cmap->element_family = {family};
  cmap->element_degree = {degree};
  cmap->create = create_{factory_name};
  cmap->geometric_dimension = {geometric_dimension};
  cmap->topological_dimension = {topological_dimension};
  cmap->is_affine = {is_affine};
  cmap->needs_permutation_data = {needs_permutation_data};
  cmap->apply_dof_transformation = apply_dof_transformation_{element_factory_name};
  cmap->get_dof_permutation = get_dof_permutation_{factory_name};
  cmap->cell_shape = {cell_shape};
  cmap->create_scalar_dofmap = create_{scalar_dofmap_name};
  return cmap;
}}

ufc_coordinate_mapping* create_coordinate_map_{prefix}(void)
{{
  return create_{factory_name}();
}}

// End of code for coordinate mapping {factory_name}
"""
