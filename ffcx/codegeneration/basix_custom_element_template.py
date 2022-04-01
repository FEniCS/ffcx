# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2018.

factory = """
// Code for custom element {factory_name}

{value_shape_init}
{wcoeffs_init}
{entity_transformations_init}
{npts_init}
{x_init}
{M_init}

ufcx_basix_custom_finite_element {factory_name} =
{{
  .cell_type = {cell_type},
  .degree = {degree},
  .value_shape = {value_shape},
  .wcoeffs_rows = {wcoeffs_rows},
  .wcoeffs_cols = {wcoeffs_cols},
  .wcoeffs = {wcoeffs},
  .entity_transformations_count = {entity_transformations_count},
  .entity_transformations_entities = {entity_transformations_entities},
  .entity_transformations_shapes = {entity_transformations_shapes},
  .entity_transformations = {entity_transformations},
  .npts = {npts},
  .x = {x},
  .M = {M},
  .map_type = {map_type},
  .discontinuous = {discontinuous},
  .highest_degree = {highest_degree},
  .highest_complete_degree = {highest_complete_degree}
}};

// End of code for custom element {factory_name}
"""
