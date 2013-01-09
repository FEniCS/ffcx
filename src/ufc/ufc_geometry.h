// This file provides utility functions for computing geometric quantities.
// This code is released into the public domain.
//
// The FEniCS Project (http://www.fenicsproject.org/) 2006-2013.

#ifndef __UFC_GEOMETRY_H
#define __UFC_GEOMETRY_H

/// Compute Jacobian of mapping from UFC reference triangle to given triangle
void compute_jacobian_triangle(double& J_00, double& J_01,
                               double& J_10, double& J_11,
                               const double* x_0,
                               const double* x_1,
                               const double* x_2)
{
  J_00 = x_1[0] - x_0[0];
  J_01 = x_2[0] - x_0[0];
  J_10 = x_1[1] - x_0[1];
  J_11 = x_2[1] - x_0[1];
}

#endif
