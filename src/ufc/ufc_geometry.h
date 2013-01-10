// This file provides utility functions for computing geometric quantities.
// This code is released into the public domain.
//
// The FEniCS Project (http://www.fenicsproject.org/) 2013.

#ifndef __UFC_GEOMETRY_H
#define __UFC_GEOMETRY_H

/// A note regarding data structures. All matrices are represented as
/// row-major flattened arrays and stored as std::vector. Benchmarks
/// indicate that with optimization (-O1 and up), std::vector is just
/// as fast as raw C++ arrays. Benchmarks also indicate that flattened
/// arrays are approximately twice as fast as nested arrays.

///--- Computation of Jacobian matrices ---

/// Compute Jacobian for 2D triangle
inline void compute_jacobian_2d(std::vector<double>& J,
                                const std::vector<double>& x)
{
  J[0] = x[2] - x[0];
  J[1] = x[4] - x[0];
  J[2] = x[3] - x[1];
  J[3] = x[5] - x[1];
}

//--- Computation of determinants ---

/// Compute determinant of 2 x 2 matrix
inline double compute_matrix_determinant_22(const std::vector<double>& A)
{
  return A[0]*A[3] - A[1]*A[2];
}

//--- Computation of matrix inverses ---

inline void compute_matrix_inverse_22(std::vector<double>& B,
                                      const std::vector<double>& A,
                                      double det)
{
  B[0] =  A[3] / det;
  B[1] = -A[1] / det;
  B[2] = -A[2] / det;
  B[3] =  A[0] / det;
}

#endif
