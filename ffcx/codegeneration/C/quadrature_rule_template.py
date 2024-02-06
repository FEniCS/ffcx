# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2023.
"""Code generation strings for a quadrature rule."""

factory = """
// Code for quadrature rule {factory_name}

{points_init}
{weights_init}

ufcx_quadrature_rule {factory_name} =
{{
  .cell_shape = {cell_shape},
  .npts = {npts},
  .topological_dimension = {topological_dimension},
  .points = {points},
  .weights = {weights}
}};

// End of code for quadrature rule {factory_name}
"""
