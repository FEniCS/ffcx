# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2018.
"""Code generation strings for a Basix custom element."""

factory = """
// Code for custom element {factory_name}

{value_shape_init}
{wcoeffs_init}
{npts_init}
{ndofs_init}
{x_init}
{M_init}

ufcx_basix_custom_finite_element {factory_name} =
{{
  .cell_type = {cell_type},
  .value_shape_length = {value_shape_length},
  .value_shape = {value_shape},
  .wcoeffs_rows = {wcoeffs_rows},
  .wcoeffs_cols = {wcoeffs_cols},
  .wcoeffs = {wcoeffs},
  .npts = {npts},
  .ndofs = {ndofs},
  .x = {x},
  .M = {M},
  .map_type = {map_type},
  .sobolev_space = {sobolev_space},
  .discontinuous = {discontinuous},
  .embedded_subdegree = {embedded_subdegree},
  .interpolation_nderivs = {interpolation_nderivs},
  .embedded_superdegree = {embedded_superdegree},
  .polyset_type = {polyset_type}
}};

// End of code for custom element {factory_name}
"""
