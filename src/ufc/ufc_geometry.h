// This file provides utility functions for computing geometric quantities.
// This code is released into the public domain.
//
// The FEniCS Project (http://www.fenicsproject.org/) 2013.

#ifndef __UFC_GEOMETRY_H
#define __UFC_GEOMETRY_H

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
/// 3. Defining an array by 'std::vector<double> x(n)', where n is a
///    literal, leads to dynamic allocation and results in significant
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
  // TODO: Move computation of det to a separate function, det is often needed when K is not
  det = J[0];
  K[0] = 1.0 / det;
}

/// Compute Jacobian (pseudo)inverse K for interval embedded in R^2
inline void compute_jacobian_inverse_interval_2d(double* K,
                                                 double& det,
                                                 const double* J)
{
  // TODO: Move computation of det to a separate function, det is often needed when K is not
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
  // TODO: Move computation of det to a separate function, det is often needed when K is not
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
  // TODO: Move computation of det to a separate function, det is often needed when K is not
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
  // TODO: Move computation of det to a separate function, det is often needed when K is not
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
  // TODO: Move computation of det to a separate function, det is often needed when K is not
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

/// Compute facet scaling factor for interval embedded in R^1
inline void compute_facet_scaling_factor_interval_1d(double & det,
                                                     const double * vertex_coordinates,
                                                     std::size_t facet)
{
  // Including this just for completeness...
  det = 1.0;
}

/// Compute facet scaling factor for interval embedded in R^2
inline void compute_facet_scaling_factor_interval_2d(double & det,
                                                     const double * vertex_coordinates,
                                                     std::size_t facet)
{
  // Including this just for completeness...
  det = 1.0;
}

/// Compute facet scaling factor for interval embedded in R^3
inline void compute_facet_scaling_factor_interval_3d(double & det,
                                                     const double * vertex_coordinates,
                                                     std::size_t facet)
{
  // Including this just for completeness...
  det = 1.0;
}

/// Compute facet scaling factor for triangle embedded in R^2
inline void compute_facet_scaling_factor_triangle_2d(double & det,
                                                     const double * vertex_coordinates,
                                                     std::size_t facet)
{
  // Get vertices on edge
  static const unsigned int edge_vertices[3][2] = {{1, 2}, {0, 2}, {0, 1}};
  const unsigned int v0 = edge_vertices[facet][0];
  const unsigned int v1 = edge_vertices[facet][1];

  // TODO: Make computation of dx* a separate function, needed for other computations as well
  // Compute scale factor (length of edge scaled by length of reference interval)
  const double dx0 = vertex_coordinates[2*v1 + 0] - vertex_coordinates[2*v0 + 0];
  const double dx1 = vertex_coordinates[2*v1 + 1] - vertex_coordinates[2*v0 + 1];

  det = std::sqrt(dx0*dx0 + dx1*dx1);
}

/// Compute facet scaling factor for triangle embedded in R^3
inline void compute_facet_scaling_factor_triangle_3d(double & det,
                                                     const double * vertex_coordinates,
                                                     std::size_t facet)
{
  // Facet determinant 2D in 3D (edge)
  // Get vertices on edge
  static const unsigned int edge_vertices[3][2] = {{1, 2}, {0, 2}, {0, 1}};
  const unsigned int v0 = edge_vertices[facet][0];
  const unsigned int v1 = edge_vertices[facet][1];

  // TODO: Make computation of dx* a separate function, needed for other computations as well
  // Compute scale factor (length of edge scaled by length of reference interval)
  const double dx0 = vertex_coordinates[3*v1 + 0] - vertex_coordinates[3*v0 + 0];
  const double dx1 = vertex_coordinates[3*v1 + 1] - vertex_coordinates[3*v0 + 1];
  const double dx2 = vertex_coordinates[3*v1 + 2] - vertex_coordinates[3*v0 + 2];

  det = std::sqrt(dx0*dx0 + dx1*dx1 + dx2*dx2);
}

/// Compute facet scaling factor for tetrahedron embedded in R^3
inline void compute_facet_scaling_factor_tetrahedron_3d(double & det,
                                                        const double * vertex_coordinates,
                                                        std::size_t facet)
{
  // Get vertices on face
  static const unsigned int face_vertices[4][3] = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};
  const unsigned int v0 = face_vertices[facet][0];
  const unsigned int v1 = face_vertices[facet][1];
  const unsigned int v2 = face_vertices[facet][2];

  // TODO: Make computation of a* a separate function, needed for other computations as well
  // Compute scale factor (area of face scaled by area of reference triangle)
  const double a0 = (vertex_coordinates[3*v0 + 1]*vertex_coordinates[3*v1 + 2]  +
                     vertex_coordinates[3*v0 + 2]*vertex_coordinates[3*v2 + 1]  +
                     vertex_coordinates[3*v1 + 1]*vertex_coordinates[3*v2 + 2]) -
                    (vertex_coordinates[3*v2 + 1]*vertex_coordinates[3*v1 + 2]  +
                     vertex_coordinates[3*v2 + 2]*vertex_coordinates[3*v0 + 1]  +
                     vertex_coordinates[3*v1 + 1]*vertex_coordinates[3*v0 + 2]);

  const double a1 = (vertex_coordinates[3*v0 + 2]*vertex_coordinates[3*v1 + 0]  +
                     vertex_coordinates[3*v0 + 0]*vertex_coordinates[3*v2 + 2]  +
                     vertex_coordinates[3*v1 + 2]*vertex_coordinates[3*v2 + 0]) -
                    (vertex_coordinates[3*v2 + 2]*vertex_coordinates[3*v1 + 0]  +
                     vertex_coordinates[3*v2 + 0]*vertex_coordinates[3*v0 + 2]  +
                     vertex_coordinates[3*v1 + 2]*vertex_coordinates[3*v0 + 0]);

  const double a2 = (vertex_coordinates[3*v0 + 0]*vertex_coordinates[3*v1 + 1]  +
                     vertex_coordinates[3*v0 + 1]*vertex_coordinates[3*v2 + 0]  +
                     vertex_coordinates[3*v1 + 0]*vertex_coordinates[3*v2 + 1]) -
                    (vertex_coordinates[3*v2 + 0]*vertex_coordinates[3*v1 + 1]  +
                     vertex_coordinates[3*v2 + 1]*vertex_coordinates[3*v0 + 0]  +
                     vertex_coordinates[3*v1 + 0]*vertex_coordinates[3*v0 + 1]);

  det = std::sqrt(a0*a0 + a1*a1 + a2*a2);
}

/// Compute facet direction for interval embedded in R^1
inline void compute_facet_normal_direction_interval_1d(bool & direction,
                                                       const double vertex_coordinates[2],
                                                       std::size_t facet)
{
  direction = facet == 0 ? vertex_coordinates[0] > vertex_coordinates[1] : vertex_coordinates[1] > vertex_coordinates[0];
}

/// Compute facet direction for triangle embedded in R^2
inline void compute_facet_normal_direction_triangle_2d(bool & direction,
                                                       const double vertex_coordinates[6],
                                                       const double dx[2],
                                                       std::size_t v0,
                                                       std::size_t facet)
{
  // FIXME: Need dx0, dx1, v0 exported from scaling factor code to use this
  direction = dx[1]*(vertex_coordinates[2*facet    ] - vertex_coordinates[2*v0    ])
            - dx[0]*(vertex_coordinates[2*facet + 1] - vertex_coordinates[2*v0 + 1])
            < 0;
}

/// Compute facet direction for tetrahedron embedded in R^3
inline void compute_facet_normal_direction_tetrahedron_3d(bool & direction,
                                                          const double vertex_coordinates[9],
                                                          const double a[3],
                                                          std::size_t v0,
                                                          std::size_t facet)
{
  // FIXME: Need a0, a1, a2, v0 exported from scaling factor code to use this
  direction = a[0]*(vertex_coordinates[3*facet    ] - vertex_coordinates[3*v0    ])
            + a[1]*(vertex_coordinates[3*facet + 1] - vertex_coordinates[3*v0 + 1])
            + a[2]*(vertex_coordinates[3*facet + 2] - vertex_coordinates[3*v0 + 2])
            < 0;
}

// TODO: Should signatures of compute_facet_normal_foo match? The snippets use different quantities.
/// Compute facet normal for interval embedded in R^1
inline void compute_facet_normal_interval_1d(double n[1],
                                             bool direction)
{
  // Facet normals are 1.0 or -1.0:   (-1.0) <-- X------X --> (1.0)
  n[0] = direction ? 1.0 : -1.0; // FIXME: Not considering facet? Got this from ffc codesnippets.
}

/// Compute facet normal for triangle embedded in R^2
inline void compute_facet_normal_triangle_2d(double n[2],
                                             const double dx[2],
                                             const double det,
                                             bool direction)
{
  // Compute facet normals from the facet scale factor constants
  n[0] = direction ?  dx[1] / det : -dx[1] / det;
  n[1] = direction ? -dx[0] / det :  dx[0] / det;
}

/// Compute facet normal for interval embedded in R^2
inline void compute_facet_normal_interval_2d(double n[2],
                                             const double vertex_coordinates[4],
                                             std::size_t facet)
{
  if (facet == 0)
  {
    n[0] = vertex_coordinates[0] - vertex_coordinates[2];
    n[1] = vertex_coordinates[1] - vertex_coordinates[3];
  }
  else
  {
    n[0] = vertex_coordinates[2] - vertex_coordinates[0];
    n[1] = vertex_coordinates[3] - vertex_coordinates[1];
  }
  const double n_length = std::sqrt(n[0]*n[0] + n[1]*n[1]);
  n[0] /= n_length;
  n[1] /= n_length;
}

/// Compute facet normal for tetrahedron embedded in R^3
inline void compute_facet_normal_tetrahedron_3d(double n[3],
                                                const double a[3],
                                                const double det,
                                                bool direction)
{
  // Compute facet normals from the facet scale factor constants
  n[0] = direction ? a[0] / det : -a[0] / det;
  n[1] = direction ? a[1] / det : -a[1] / det;
  n[2] = direction ? a[2] / det : -a[2] / det;
}

/// Compute facet normal for triangle embedded in R^3
inline void compute_facet_normal_triangle_3d(double n[3],
                                             const double vertex_coordinates[6],
                                             std::size_t facet)
{
  // Compute facet normal for triangles in 3D
  const unsigned int vertex0 = facet;

  // Get coordinates corresponding the vertex opposite this
  static const unsigned int edge_vertices[3][2] = {{1, 2}, {0, 2}, {0, 1}};
  const unsigned int vertex1 = edge_vertices[facet][0];
  const unsigned int vertex2 = edge_vertices[facet][1];

  // Define vectors n = (p2 - p0) and t = normalized (p2 - p1)
  n[0] = vertex_coordinates[3*vertex2 + 0] - vertex_coordinates[3*vertex0 + 0];
  n[1] = vertex_coordinates[3*vertex2 + 1] - vertex_coordinates[3*vertex0 + 1];
  n[2] = vertex_coordinates[3*vertex2 + 2] - vertex_coordinates[3*vertex0 + 2];

  double t0 = vertex_coordinates[3*vertex2 + 0] - vertex_coordinates[3*vertex1 + 0];
  double t1 = vertex_coordinates[3*vertex2 + 1] - vertex_coordinates[3*vertex1 + 1];
  double t2 = vertex_coordinates[3*vertex2 + 2] - vertex_coordinates[3*vertex1 + 2];
  const double t_length = std::sqrt(t0*t0 + t1*t1 + t2*t2);
  t0 /= t_length;
  t1 /= t_length;
  t2 /= t_length;

  // Subtract, the projection of (p2  - p0) onto (p2 - p1), from (p2 - p0)
  const double ndott = t0*n[0] + t1*n[1] + t2*n[2];
  n[0] -= ndott*t0;
  n[1] -= ndott*t1;
  n[2] -= ndott*t2;
  const double n_length = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);

  // Normalize
  n[0] /= n_length;
  n[1] /= n_length;
  n[2] /= n_length;
}

/// Compute facet normal for interval embedded in R^3
inline void compute_facet_normal_interval_3d(double n[3],
                                             const double vertex_coordinates[6],
                                             std::size_t facet)
{
  if (facet == 0)
  {
    n[0] = vertex_coordinates[0] - vertex_coordinates[3];
    n[1] = vertex_coordinates[1] - vertex_coordinates[4];
    n[1] = vertex_coordinates[2] - vertex_coordinates[5];
  }
  else
  {
    n[0] = vertex_coordinates[3] - vertex_coordinates[0];
    n[1] = vertex_coordinates[4] - vertex_coordinates[1];
    n[1] = vertex_coordinates[5] - vertex_coordinates[2];
  }
  const double n_length = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  n[0] /= n_length;
  n[1] /= n_length;
  n[2] /= n_length;
}

#endif
