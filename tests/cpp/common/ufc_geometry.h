// This file provides utility functions for computing geometric quantities.
// This code is released into the public domain.
//
// The FEniCS Project (http://www.fenicsproject.org/) 2013.

#ifndef __UFLACS_GEOMETRY_H
#define __UFLACS_GEOMETRY_H

#include <cmath>

/// A note regarding data structures. All matrices are represented as
/// row-major flattened raw C++ arrays. Benchmarks indicate that when
/// optimization (-O1 and up) is used, the following conditions hold:
///
/// 1. std::vector is just as fast as raw C++ arrays for indexing.
///
/// 2. Flattened arrays are twice as fast as nested arrays, both for
///    std:vector and raw C++ arrays.
///
/// 3. Defining an array by 'std::vector<double> x(n)' where n is a
///    literal leads to dynamic allocatioin and results in significant
///    slowdowns compared to the definition 'double x[n]'.
///
/// The conclusion is that we should use flattened raw C++ arrays in
/// the interfaces for these utility functions, since some of the
/// arrays passed to these functions (in particular Jacobians) are
/// created inside the generated functions (tabulate_tensor). Note
/// that an std::vector x may also be passed as raw pointer by &x[0].

///--- Computation of Jacobian matrices ---

/// Compute Jacobian J for interval embedded in R^1
inline void compute_jacobian_interval_1d(double* J,
                                         const double* vertex_coordinates)
{
  J[0] = vertex_coordinates[1] - vertex_coordinates[0];
}

/// Compute Jacobian J for interval embedded in R^2
inline void compute_jacobian_interval_2d(double* J,
                                         const double* vertex_coordinates)
{
  J[0] = vertex_coordinates[2] - vertex_coordinates[0];
  J[1] = vertex_coordinates[3] - vertex_coordinates[1];
}

/// Compute Jacobian J for interval embedded in R^3
inline void compute_jacobian_interval_3d(double* J,
                                         const double* vertex_coordinates)
{
  J[0] = vertex_coordinates[3] - vertex_coordinates[0];
  J[1] = vertex_coordinates[4] - vertex_coordinates[1];
  J[2] = vertex_coordinates[5] - vertex_coordinates[2];
}

/// Compute Jacobian J for triangle embedded in R^2
inline void compute_jacobian_triangle_2d(double* J,
                                         const double* vertex_coordinates)
{
  J[0] = vertex_coordinates[2] - vertex_coordinates[0];
  J[1] = vertex_coordinates[4] - vertex_coordinates[0];
  J[2] = vertex_coordinates[3] - vertex_coordinates[1];
  J[3] = vertex_coordinates[5] - vertex_coordinates[1];
}

/// Compute Jacobian J for triangle embedded in R^3
inline void compute_jacobian_triangle_3d(double* J,
                                         const double* vertex_coordinates)
{
  J[0] = vertex_coordinates[3] - vertex_coordinates[0];
  J[1] = vertex_coordinates[6] - vertex_coordinates[0];
  J[2] = vertex_coordinates[4] - vertex_coordinates[1];
  J[3] = vertex_coordinates[7] - vertex_coordinates[1];
  J[4] = vertex_coordinates[5] - vertex_coordinates[2];
  J[5] = vertex_coordinates[8] - vertex_coordinates[2];
}

/// Compute Jacobian J for tetrahedron embedded in R^3
inline void compute_jacobian_tetrahedron_3d(double* J,
                                            const double* vertex_coordinates)
{
  J[0] = vertex_coordinates[3]  - vertex_coordinates[0];
  J[1] = vertex_coordinates[6]  - vertex_coordinates[0];
  J[2] = vertex_coordinates[9]  - vertex_coordinates[0];
  J[3] = vertex_coordinates[4]  - vertex_coordinates[1];
  J[4] = vertex_coordinates[7]  - vertex_coordinates[1];
  J[5] = vertex_coordinates[10] - vertex_coordinates[1];
  J[6] = vertex_coordinates[5]  - vertex_coordinates[2];
  J[7] = vertex_coordinates[8]  - vertex_coordinates[2];
  J[8] = vertex_coordinates[11] - vertex_coordinates[2];
}

//--- Computation of Jacobian inverses ---

/// Compute Jacobian inverse K for interval embedded in R^1
inline void compute_jacobian_inverse_interval_1d(double* K,
                                                 double& det,
                                                 const double* J)
{
  det = J[0];
  K[0] = 1.0 / det;
}

/// Compute Jacobian (pseudo)inverse K for interval embedded in R^2
inline void compute_jacobian_inverse_interval_2d(double* K,
                                                 double& det,
                                                 const double* J)
{
  const double det2 = J[0]*J[0] + J[1]*J[1];
  det = std::sqrt(det2);

  K[0] = J[0] / det2;
  K[1] = J[1] / det2;
}

/// Compute Jacobian (pseudo)inverse K for interval embedded in R^3
inline void compute_jacobian_inverse_interval_3d(double* K,
                                                 double& det,
                                                 const double* J)
{
  const double det2 = J[0]*J[0] + J[1]*J[1] + J[2]*J[2];
  det = std::sqrt(det2);

  K[0] = J[0] / det2;
  K[1] = J[1] / det2;
  K[2] = J[2] / det2;
}

/// Compute Jacobian inverse K for triangle embedded in R^2
inline void compute_jacobian_inverse_triangle_2d(double* K,
                                                 double& det,
                                                 const double* J)
{
  det = J[0]*J[3] - J[1]*J[2];

  K[0] =  J[3] / det;
  K[1] = -J[1] / det;
  K[2] = -J[2] / det;
  K[3] =  J[0] / det;
}

/// Compute Jacobian (pseudo)inverse K for triangle embedded in R^3
inline void compute_jacobian_inverse_triangle_3d(double* K,
                                                 double& det,
                                                 const double* J)
{
  const double d_0 = J[2]*J[5] - J[4]*J[3];
  const double d_1 = J[4]*J[1] - J[0]*J[5];
  const double d_2 = J[0]*J[3] - J[2]*J[1];

  const double c_0 = J[0]*J[0] + J[2]*J[2] + J[4]*J[4];
  const double c_1 = J[1]*J[1] + J[3]*J[3] + J[5]*J[5];
  const double c_2 = J[0]*J[1] + J[2]*J[3] + J[4]*J[5];

  const double den = c_0*c_1 - c_2*c_2;

  const double det2 = d_0*d_0 + d_1*d_1 + d_2*d_2;
  det = std::sqrt(det2);

  K[0] = (J[0]*c_1 - J[1]*c_2) / den;
  K[1] = (J[2]*c_1 - J[3]*c_2) / den;
  K[2] = (J[4]*c_1 - J[5]*c_2) / den;
  K[3] = (J[1]*c_0 - J[0]*c_2) / den;
  K[4] = (J[3]*c_0 - J[2]*c_2) / den;
  K[5] = (J[5]*c_0 - J[4]*c_2) / den;
}

/// Compute Jacobian inverse K for tetrahedron embedded in R^3
inline void compute_jacobian_inverse_tetrahedron_3d(double* K,
                                                    double& det,
                                                    const double* J)
{
  const double d_00 = J[4]*J[8] - J[5]*J[7];
  const double d_01 = J[5]*J[6] - J[3]*J[8];
  const double d_02 = J[3]*J[7] - J[4]*J[6];
  const double d_10 = J[2]*J[7] - J[1]*J[8];
  const double d_11 = J[0]*J[8] - J[2]*J[6];
  const double d_12 = J[1]*J[6] - J[0]*J[7];
  const double d_20 = J[1]*J[5] - J[2]*J[4];
  const double d_21 = J[2]*J[3] - J[0]*J[5];
  const double d_22 = J[0]*J[4] - J[1]*J[3];

  det = J[0]*d_00 + J[3]*d_10 + J[6]*d_20;

  K[0] = d_00 / det;
  K[1] = d_10 / det;
  K[2] = d_20 / det;
  K[3] = d_01 / det;
  K[4] = d_11 / det;
  K[5] = d_21 / det;
  K[6] = d_02 / det;
  K[7] = d_12 / det;
  K[8] = d_22 / det;
}

#endif
