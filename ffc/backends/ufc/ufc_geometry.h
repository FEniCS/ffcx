// This file provides utility functions for computing geometric quantities.
// This code is released into the public domain.
//
// The FEniCS Project (http://www.fenicsproject.org/) 2013-2015.

#ifndef __UFC_GEOMETRY_H
#define __UFC_GEOMETRY_H

#include <cmath>

// TODO: Wrap all in namespace ufc
//namespace ufc
//{

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

// TODO: Should signatures of compute_<foo>_<cell>_<n>d match for each foo?
//       On one hand the snippets use different quantities, on the other
//       some consistency is nice to simplify the code generation.
//       Currently only the arguments that are actually used are included.

/// --- Some fixed numbers by name for readability in this file ---
// TODO: Use these numbers where relevant below to make this file more self documenting

// (namespaced using UFC_ in the names because they collide with variables in other libraries)

// Use this for array dimensions indexed by local vertex number
#define UFC_NUM_VERTICES_IN_INTERVAL 2
#define UFC_NUM_VERTICES_IN_TRIANGLE 3
#define UFC_NUM_VERTICES_IN_TETRAHEDRON 4
#define UFC_NUM_VERTICES_IN_QUADRILATERAL 4
#define UFC_NUM_VERTICES_IN_HEXAHEDRON 8

// Use this for array dimensions indexed by local edge number
#define UFC_NUM_EDGES_IN_INTERVAL 1
#define UFC_NUM_EDGES_IN_TRIANGLE 3
#define UFC_NUM_EDGES_IN_TETRAHEDRON 6
#define UFC_NUM_EDGES_IN_QUADRILATERAL 4
#define UFC_NUM_EDGES_IN_HEXAHEDRON 12

// Use this for array dimensions indexed by local facet number
#define UFC_NUM_FACETS_IN_INTERVAL 2
#define UFC_NUM_FACETS_IN_TRIANGLE 3
#define UFC_NUM_FACETS_IN_TETRAHEDRON 4
#define UFC_NUM_FACETS_IN_QUADRILATERAL 4
#define UFC_NUM_FACETS_IN_HEXAHEDRON 6

// Use UFC_GDIM_N to show the intention that the geometric dimension is N
#define UFC_GDIM_1 1
#define UFC_GDIM_2 2
#define UFC_GDIM_3 3

// Use UFC_TDIM_N to show the intention that the topological dimension is N
#define UFC_TDIM_1 1
#define UFC_TDIM_2 2
#define UFC_TDIM_3 3


/// --- Local reference cell coordinates by UFC conventions ---

static const double interval_vertices[UFC_NUM_VERTICES_IN_INTERVAL][UFC_TDIM_1] = {
  {0.0},
  {1.0}
  };

static const double triangle_vertices[UFC_NUM_VERTICES_IN_TRIANGLE][UFC_TDIM_2] = {
  {0.0, 0.0},
  {1.0, 0.0},
  {0.0, 1.0}
  };

static const double tetrahedron_vertices[UFC_NUM_VERTICES_IN_TETRAHEDRON][UFC_TDIM_3] = {
  {0.0, 0.0, 0.0},
  {1.0, 0.0, 0.0},
  {0.0, 1.0, 0.0},
  {0.0, 0.0, 1.0}
  };

static const double quadrilateral_vertices[UFC_NUM_VERTICES_IN_QUADRILATERAL][UFC_TDIM_2] = {
  {0.0, 0.0},
  {0.0, 1.0},
  {1.0, 0.0},
  {1.0, 1.0},
  };

static const double hexahedron_vertices[UFC_NUM_VERTICES_IN_HEXAHEDRON][UFC_TDIM_3] = {
  {0.0, 0.0, 0.0},
  {0.0, 0.0, 1.0},
  {0.0, 1.0, 0.0},
  {0.0, 1.0, 1.0},
  {1.0, 0.0, 0.0},
  {1.0, 0.0, 1.0},
  {1.0, 1.0, 0.0},
  {1.0, 1.0, 1.0},
  };

/// --- Local reference cell midpoint by UFC conventions ---

static const double interval_midpoint[UFC_TDIM_1] = {
  0.5
  };

static const double triangle_midpoint[UFC_TDIM_2] = {
  1.0/3.0, 1.0/3.0
  };

static const double tetrahedron_midpoint[UFC_TDIM_3] = {
  0.25, 0.25, 0.25
  };

static const double quadrilateral_midpoint[UFC_TDIM_2] = {
  0.5, 0.5
  };

static const double hexahedron_midpoint[UFC_TDIM_3] = {
  0.5, 0.5, 0.5
  };

/// --- Local reference cell facet midpoints by UFC conventions ---

static const double interval_facet_midpoint[UFC_NUM_FACETS_IN_INTERVAL][UFC_TDIM_1] = {
  {0.0},
  {1.0}
  };

static const double triangle_facet_midpoint[UFC_NUM_FACETS_IN_TRIANGLE][UFC_TDIM_2] = {
  {0.5, 0.5},
  {0.0, 0.5},
  {0.5, 0.0}
  };

static const double tetrahedron_facet_midpoint[UFC_NUM_FACETS_IN_TETRAHEDRON][UFC_TDIM_3] = {
  {0.5, 0.5, 0.5},
  {0.0, 1.0/3.0, 1.0/3.0},
  {1.0/3.0, 0.0, 1.0/3.0},
  {1.0/3.0, 1.0/3.0, 0.0},
  };

static const double quadrilateral_facet_midpoint[UFC_NUM_FACETS_IN_QUADRILATERAL][UFC_TDIM_2] = {
  {0.0, 0.5},
  {1.0, 0.5},
  {0.5, 0.0},
  {0.5, 1.0},
  };

static const double hexahedron_facet_midpoint[UFC_NUM_FACETS_IN_HEXAHEDRON][UFC_TDIM_3] = {
  {0.0, 0.5, 0.5},
  {1.0, 0.5, 0.5},
  {0.5, 0.0, 0.5},
  {0.5, 1.0, 0.5},
  {0.5, 0.5, 0.0},
  {0.5, 0.5, 1.0},
  };

/// --- Local reference cell facet orientations by UFC conventions ---

static const double interval_facet_orientations[UFC_NUM_FACETS_IN_INTERVAL] = {
  -1.0,
  +1.0,
  };

static const double triangle_facet_orientations[UFC_NUM_FACETS_IN_TRIANGLE] = {
  +1.0,
  -1.0,
  +1.0
  };

static const double tetrahedron_facet_orientations[UFC_NUM_FACETS_IN_TETRAHEDRON] = {
  +1.0,
  -1.0,
  +1.0,
  -1.0
  };

// FIXME: Insert quad conventions here
/*
static const double quadrilateral_facet_orientations[UFC_NUM_FACETS_IN_QUADRILATERAL] = {
  +1.0,
  -1.0,
  -1.0,
  +1.0,
  };
*/

// FIXME: Insert quad conventions here
/*
static const double hexahedron_facet_orientations[UFC_NUM_FACETS_IN_HEXAHEDRON] = {
  +1.0,
  -1.0,
  -1.0,
  +1.0,
  +1.0,
  -1.0,
  };
*/

/// --- Local reference cell entity relations by UFC conventions ---

static const unsigned int triangle_edge_vertices[UFC_NUM_EDGES_IN_TRIANGLE][2] = {
  {1, 2},
  {0, 2},
  {0, 1}
  };

static const unsigned int tetrahedron_edge_vertices[UFC_NUM_EDGES_IN_TETRAHEDRON][2] = {
  {2, 3},
  {1, 3},
  {1, 2},
  {0, 3},
  {0, 2},
  {0, 1}
  };

static const unsigned int quadrilateral_edge_vertices[UFC_NUM_EDGES_IN_QUADRILATERAL][2] = {
  {0, 1},
  {2, 3},
  {0, 2},
  {1, 3},
  };

static const unsigned int hexahedron_edge_vertices[UFC_NUM_EDGES_IN_HEXAHEDRON][2] = {
  {0, 1},
  {2, 3},
  {4, 5},
  {6, 7},
  {0, 2},
  {1, 3},
  {4, 6},
  {5, 7},
  {0, 4},
  {1, 5},
  {2, 6},
  {3, 7},
  };


/// --- Local reference cell entity relations by UFC conventions ---

static const unsigned int interval_facet_vertices[UFC_NUM_FACETS_IN_INTERVAL][1] = {
  {0},
  {1}
  };

static const unsigned int triangle_facet_vertices[UFC_NUM_FACETS_IN_TRIANGLE][UFC_NUM_VERTICES_IN_INTERVAL] = {
  {1, 2},
  {0, 2},
  {0, 1}
  };

static const unsigned int tetrahedron_facet_vertices[UFC_NUM_FACETS_IN_TETRAHEDRON][UFC_NUM_VERTICES_IN_TRIANGLE] = {
  {1, 2, 3},
  {0, 2, 3},
  {0, 1, 3},
  {0, 1, 2}
  };

static const unsigned int tetrahedron_facet_edge_vertices[UFC_NUM_FACETS_IN_TETRAHEDRON][UFC_NUM_FACETS_IN_TRIANGLE][UFC_NUM_VERTICES_IN_INTERVAL] = {
  {{2, 3}, {1, 3}, {1, 2}},
  {{2, 3}, {0, 3}, {0, 2}},
  {{1, 3}, {0, 3}, {0, 1}},
  {{1, 2}, {0, 2}, {0, 1}}
  };

static const unsigned int quadrilateral_facet_vertices[UFC_NUM_FACETS_IN_QUADRILATERAL][UFC_NUM_VERTICES_IN_INTERVAL] = {
  {0, 1},
  {2, 3},
  {0, 2},
  {1, 3},
  };

static const unsigned int hexahedron_facet_vertices[UFC_NUM_FACETS_IN_HEXAHEDRON][UFC_NUM_VERTICES_IN_QUADRILATERAL] = {
  {0, 1, 2, 3},
  {4, 5, 6, 7},
  {0, 1, 4, 5},
  {2, 3, 6, 7},
  {0, 2, 4, 6},
  {1, 3, 5, 7},
  };

static const unsigned int hexahedron_facet_edge_vertices[UFC_NUM_FACETS_IN_HEXAHEDRON][UFC_NUM_FACETS_IN_QUADRILATERAL][UFC_NUM_VERTICES_IN_INTERVAL] = {
  {{0, 1}, {2, 3}, {0, 2}, {1, 3}},
  {{4, 5}, {6, 7}, {4, 6}, {5, 7}},
  {{0, 1}, {4, 5}, {0, 4}, {1, 5}},
  {{2, 3}, {6, 7}, {2, 6}, {3, 7}},
  {{0, 2}, {4, 6}, {0, 4}, {2, 6}},
  {{1, 3}, {5, 7}, {1, 5}, {3, 7}},
  };

/// --- Reference cell edge vectors by UFC conventions (edge vertex 1 - edge vertex 0 for each edge in cell) ---

static const double triangle_reference_edge_vectors[UFC_NUM_EDGES_IN_TRIANGLE][UFC_TDIM_2] = {
  {-1.0, 1.0},
  { 0.0, 1.0},
  { 1.0, 0.0},
  };

static const double tetrahedron_reference_edge_vectors[UFC_NUM_EDGES_IN_TETRAHEDRON][UFC_TDIM_3] = {
  { 0.0, -1.0,  1.0},
  {-1.0,  0.0,  1.0},
  {-1.0,  1.0,  0.0},
  { 0.0,  0.0,  1.0},
  { 0.0,  1.0,  0.0},
  { 1.0,  0.0,  0.0},
  };

// Edge vectors for each triangle facet of a tetrahedron
static const double tetrahedron_facet_reference_edge_vectors[UFC_NUM_FACETS_IN_TETRAHEDRON][UFC_NUM_EDGES_IN_TRIANGLE][UFC_TDIM_3] = {
  { // facet 0
    { 0.0, -1.0,  1.0},
    {-1.0,  0.0,  1.0},
    {-1.0,  1.0,  0.0},
  },
  { // facet 1
    { 0.0, -1.0,  1.0},
    { 0.0,  0.0,  1.0},
    { 0.0,  1.0,  0.0},
  },
  { // facet 2
    {-1.0,  0.0,  1.0},
    { 0.0,  0.0,  1.0},
    { 1.0,  0.0,  0.0},
  },
  { // facet 3
    {-1.0,  1.0,  0.0},
    { 0.0,  1.0,  0.0},
    { 1.0,  0.0,  0.0},
  },
  };

static const double quadrilateral_reference_edge_vectors[UFC_NUM_EDGES_IN_QUADRILATERAL][UFC_TDIM_2] = {
  { 0.0, 1.0},
  { 0.0, 1.0},
  { 1.0, 0.0},
  { 1.0, 0.0},
  };

static const double hexahedron_reference_edge_vectors[UFC_NUM_EDGES_IN_HEXAHEDRON][UFC_TDIM_3] = {
  { 0.0,  0.0,  1.0},
  { 0.0,  0.0,  1.0},
  { 0.0,  0.0,  1.0},
  { 0.0,  0.0,  1.0},
  { 0.0,  1.0,  0.0},
  { 0.0,  1.0,  0.0},
  { 0.0,  1.0,  0.0},
  { 0.0,  1.0,  0.0},
  { 1.0,  0.0,  0.0},
  { 1.0,  0.0,  0.0},
  { 1.0,  0.0,  0.0},
  { 1.0,  0.0,  0.0},
  };

// Edge vectors for each quadrilateral facet of a hexahedron
static const double hexahedron_facet_reference_edge_vectors[UFC_NUM_FACETS_IN_HEXAHEDRON][UFC_NUM_EDGES_IN_QUADRILATERAL][UFC_TDIM_3] = {
  { // facet 0
    { 0.0,  0.0,  1.0},
    { 0.0,  0.0,  1.0},
    { 0.0,  1.0,  0.0},
    { 0.0,  1.0,  0.0},
  },
  { // facet 1
    { 0.0,  0.0,  1.0},
    { 0.0,  0.0,  1.0},
    { 0.0,  1.0,  0.0},
    { 0.0,  1.0,  0.0},
  },
  { // facet 2
    { 0.0,  0.0,  1.0},
    { 0.0,  0.0,  1.0},
    { 1.0,  0.0,  0.0},
    { 1.0,  0.0,  0.0},
  },
  { // facet 3
    { 0.0,  0.0,  1.0},
    { 0.0,  0.0,  1.0},
    { 1.0,  0.0,  0.0},
    { 1.0,  0.0,  0.0},
  },
  { // facet 4
    { 0.0,  1.0,  0.0},
    { 0.0,  1.0,  0.0},
    { 1.0,  0.0,  0.0},
    { 1.0,  0.0,  0.0},
  },
  { // facet 5
    { 0.0,  1.0,  0.0},
    { 0.0,  1.0,  0.0},
    { 1.0,  0.0,  0.0},
    { 1.0,  0.0,  0.0},
  },
  };

/// --- Reference cell facet normals by UFC conventions (outwards pointing on reference cell) ---

static const double interval_reference_facet_normals[UFC_NUM_FACETS_IN_INTERVAL][UFC_TDIM_1] = {
  {-1.0},
  {+1.0},
  };

static const double triangle_reference_facet_normals[UFC_NUM_FACETS_IN_TRIANGLE][UFC_TDIM_2] = {
  { 0.7071067811865476, 0.7071067811865476 },
  {-1.0,  0.0},
  { 0.0, -1.0},
  };

static const double tetrahedron_reference_facet_normals[UFC_NUM_FACETS_IN_TETRAHEDRON][UFC_TDIM_3] = {
  {0.5773502691896258, 0.5773502691896258, 0.5773502691896258},
  {-1.0,  0.0,  0.0},
  { 0.0, -1.0,  0.0},
  { 0.0,  0.0, -1.0},
  };

static const double quadrilateral_reference_facet_normals[UFC_NUM_FACETS_IN_QUADRILATERAL][UFC_TDIM_2] = {
  { -1.0,  0.0 },
  {  1.0,  0.0 },
  {  0.0, -1.0 },
  {  0.0,  1.0 },
  };

static const double hexahedron_reference_facet_normals[UFC_NUM_FACETS_IN_HEXAHEDRON][UFC_TDIM_3] = {
  { -1.0,  0.0,  0.0},
  {  1.0,  0.0,  0.0},
  {  0.0, -1.0,  0.0},
  {  0.0,  1.0,  0.0},
  {  0.0,  0.0, -1.0},
  {  0.0,  0.0,  1.0},
  };

/// --- Reference cell volumes by UFC conventions ---

static const double interval_reference_cell_volume = 1.0;
static const double triangle_reference_cell_volume = 0.5;
static const double tetrahedron_reference_cell_volume = 1.0/6.0;
static const double quadrilateral_reference_cell_volume = 1.0;
static const double hexahedron_reference_cell_volume = 1.0;

static const double interval_reference_facet_volume = 1.0;
static const double triangle_reference_facet_volume = 1.0;
static const double tetrahedron_reference_facet_volume = 0.5;
static const double quadrilateral_reference_facet_volume = 1.0;
static const double hexahedron_reference_facet_volume = 1.0;

/// --- Jacobians of reference facet cell to reference cell coordinate mappings by UFC conventions ---

static const double triangle_reference_facet_jacobian[UFC_NUM_FACETS_IN_TRIANGLE][UFC_TDIM_2][UFC_TDIM_2-1] = {
  { {-1.0}, { 1.0} },
  { { 0.0}, { 1.0} },
  { { 1.0}, { 0.0} },
  };

static const double tetrahedron_reference_facet_jacobian[UFC_NUM_FACETS_IN_TETRAHEDRON][UFC_TDIM_3][UFC_TDIM_3-1] = {
  { {-1.0, -1.0}, {1.0, 0.0}, {0.0, 1.0} },
  { { 0.0,  0.0}, {1.0, 0.0}, {0.0, 1.0} },
  { { 1.0,  0.0}, {0.0, 0.0}, {0.0, 1.0} },
  { { 1.0,  0.0}, {0.0, 1.0}, {0.0, 0.0} },
  };

static const double quadrilateral_reference_facet_jacobian[UFC_NUM_FACETS_IN_QUADRILATERAL][UFC_TDIM_2][UFC_TDIM_2-1] = {
  { { 0.0}, { 1.0} },
  { { 0.0}, { 1.0} },
  { { 1.0}, { 0.0} },
  { { 1.0}, { 0.0} },
  };

static const double hexahedron_reference_facet_jacobian[UFC_NUM_FACETS_IN_HEXAHEDRON][UFC_TDIM_3][UFC_TDIM_3-1] = {
  { { 0.0,  0.0}, {1.0, 0.0}, {0.0, 1.0} },
  { { 0.0,  0.0}, {1.0, 0.0}, {0.0, 1.0} },
  { { 1.0,  0.0}, {0.0, 0.0}, {0.0, 1.0} },
  { { 1.0,  0.0}, {0.0, 0.0}, {0.0, 1.0} },
  { { 1.0,  0.0}, {0.0, 1.0}, {0.0, 0.0} },
  { { 1.0,  0.0}, {0.0, 1.0}, {0.0, 0.0} },
  };

/// --- Coordinate mappings from reference facet cell to reference cell by UFC conventions ---

inline void compute_reference_facet_to_reference_cell_coordinates_interval(double Xc[UFC_TDIM_1], unsigned int facet)
{
  switch (facet)
  {
  case 0:
    Xc[0] = 0.0;
    break;
  case 1:
    Xc[0] = 1.0;
    break;
  };
}

inline void compute_reference_facet_to_reference_cell_coordinates_triangle(double Xc[UFC_TDIM_2], const double Xf[UFC_TDIM_2-1], unsigned int facet)
{
  switch (facet)
  {
  case 0:
    Xc[0] = 1.0 - Xf[0];
    Xc[1] = Xf[0];
    break;
  case 1:
    Xc[0] = 0.0;
    Xc[1] = Xf[0];
    break;
  case 2:
    Xc[0] = Xf[0];
    Xc[1] = 0.0;
    break;
  };
}

inline void compute_reference_facet_to_reference_cell_coordinates_tetrahedron(double Xc[UFC_TDIM_3], const double Xf[UFC_TDIM_3-1], unsigned int facet)
{
  switch (facet)
  {
  case 0:
    Xc[0] = 1.0 - Xf[0] - Xf[1];
    Xc[1] = Xf[0];
    Xc[2] = Xf[1];
    break;
  case 1:
    Xc[0] = 0.0;
    Xc[1] = Xf[0];
    Xc[2] = Xf[1];
    break;
  case 2:
    Xc[0] = Xf[0];
    Xc[1] = 0.0;
    Xc[2] = Xf[1];
    break;
  case 3:
    Xc[0] = Xf[0];
    Xc[1] = Xf[1];
    Xc[2] = 0.0;
    break;
  };
}


///--- Computation of Jacobian matrices ---

/// Compute Jacobian J for interval embedded in R^1
inline void compute_jacobian_interval_1d(double J[UFC_GDIM_1*UFC_TDIM_1],
                                         const double coordinate_dofs[2])
{
  J[0] = coordinate_dofs[1] - coordinate_dofs[0];
}

/// Compute Jacobian J for interval embedded in R^2
inline void compute_jacobian_interval_2d(double J[UFC_GDIM_2*UFC_TDIM_1],
                                         const double coordinate_dofs[4])
{
  J[0] = coordinate_dofs[2] - coordinate_dofs[0];
  J[1] = coordinate_dofs[3] - coordinate_dofs[1];
}

/// Compute Jacobian J for interval embedded in R^3
inline void compute_jacobian_interval_3d(double J[UFC_GDIM_3*UFC_TDIM_1],
                                         const double coordinate_dofs[6])
{
  J[0] = coordinate_dofs[3] - coordinate_dofs[0];
  J[1] = coordinate_dofs[4] - coordinate_dofs[1];
  J[2] = coordinate_dofs[5] - coordinate_dofs[2];
}

/// Compute Jacobian J for triangle embedded in R^2
inline void compute_jacobian_triangle_2d(double J[UFC_GDIM_2*UFC_TDIM_2],
                                         const double coordinate_dofs[6])
{
  J[0] = coordinate_dofs[2] - coordinate_dofs[0];
  J[1] = coordinate_dofs[4] - coordinate_dofs[0];
  J[2] = coordinate_dofs[3] - coordinate_dofs[1];
  J[3] = coordinate_dofs[5] - coordinate_dofs[1];
}

/// Compute Jacobian J for triangle embedded in R^3
inline void compute_jacobian_triangle_3d(double J[UFC_GDIM_3*UFC_TDIM_2],
                                         const double coordinate_dofs[9])
{
  J[0] = coordinate_dofs[3] - coordinate_dofs[0];
  J[1] = coordinate_dofs[6] - coordinate_dofs[0];
  J[2] = coordinate_dofs[4] - coordinate_dofs[1];
  J[3] = coordinate_dofs[7] - coordinate_dofs[1];
  J[4] = coordinate_dofs[5] - coordinate_dofs[2];
  J[5] = coordinate_dofs[8] - coordinate_dofs[2];
}

/// Compute Jacobian J for tetrahedron embedded in R^3
inline void compute_jacobian_tetrahedron_3d(double J[UFC_GDIM_3*UFC_TDIM_3],
                                            const double coordinate_dofs[12])
{
  J[0] = coordinate_dofs[3]  - coordinate_dofs[0];
  J[1] = coordinate_dofs[6]  - coordinate_dofs[0];
  J[2] = coordinate_dofs[9]  - coordinate_dofs[0];
  J[3] = coordinate_dofs[4]  - coordinate_dofs[1];
  J[4] = coordinate_dofs[7]  - coordinate_dofs[1];
  J[5] = coordinate_dofs[10] - coordinate_dofs[1];
  J[6] = coordinate_dofs[5]  - coordinate_dofs[2];
  J[7] = coordinate_dofs[8]  - coordinate_dofs[2];
  J[8] = coordinate_dofs[11] - coordinate_dofs[2];
}

//--- Computation of Jacobian inverses --- // TODO: Remove this when ffc is updated to use the NEW ones below

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

//--- NEW Computation of Jacobian (sub)determinants ---

/// Compute Jacobian determinant for interval embedded in R^1
inline void compute_jacobian_determinants_interval_1d(double & det,
                                                      const double J[UFC_GDIM_1*UFC_TDIM_1])
{
  det = J[0];
}

/// Compute Jacobian (pseudo)determinants for interval embedded in R^2
inline void compute_jacobian_determinants_interval_2d(double & det2,
                                                      double & det,
                                                      const double J[UFC_GDIM_2*UFC_TDIM_1])
{
  det2 = J[0]*J[0] + J[1]*J[1];
  det = std::sqrt(det2);
}

/// Compute Jacobian (pseudo)determinants for interval embedded in R^3
inline void compute_jacobian_determinants_interval_3d(double & det2,
                                                      double & det,
                                                      const double J[UFC_GDIM_3*UFC_TDIM_1])
{
  det2 = J[0]*J[0] + J[1]*J[1] + J[2]*J[2];
  det = std::sqrt(det2);
}

/// Compute Jacobian determinant for triangle embedded in R^2
inline void compute_jacobian_determinants_triangle_2d(double & det,
                                                      const double J[UFC_GDIM_2*UFC_TDIM_2])
{
  det = J[0]*J[3] - J[1]*J[2];
}

/// Compute Jacobian (pseudo)determinants for triangle embedded in R^3
inline void compute_jacobian_determinants_triangle_3d(double & den,
                                                      double & det2,
                                                      double & det,
                                                      double c[3],
                                                      const double J[UFC_GDIM_3*UFC_TDIM_2])
{
  const double d_0 = J[2]*J[5] - J[4]*J[3];
  const double d_1 = J[4]*J[1] - J[0]*J[5];
  const double d_2 = J[0]*J[3] - J[2]*J[1];

  c[0] = J[0]*J[0] + J[2]*J[2] + J[4]*J[4];
  c[1] = J[1]*J[1] + J[3]*J[3] + J[5]*J[5];
  c[2] = J[0]*J[1] + J[2]*J[3] + J[4]*J[5];

  den = c[0]*c[1] - c[2]*c[2];

  det2 = d_0*d_0 + d_1*d_1 + d_2*d_2;
  det = std::sqrt(det2);
}

/// Compute Jacobian determinants for tetrahedron embedded in R^3
inline void compute_jacobian_determinants_tetrahedron_3d(double & det,
                                                         double d[9],
                                                         const double J[UFC_GDIM_3*UFC_TDIM_3])
{
  d[0*3 + 0] = J[4]*J[8] - J[5]*J[7];
  d[0*3 + 1] = J[5]*J[6] - J[3]*J[8];
  d[0*3 + 2] = J[3]*J[7] - J[4]*J[6];
  d[1*3 + 0] = J[2]*J[7] - J[1]*J[8];
  d[1*3 + 1] = J[0]*J[8] - J[2]*J[6];
  d[1*3 + 2] = J[1]*J[6] - J[0]*J[7];
  d[2*3 + 0] = J[1]*J[5] - J[2]*J[4];
  d[2*3 + 1] = J[2]*J[3] - J[0]*J[5];
  d[2*3 + 2] = J[0]*J[4] - J[1]*J[3];

  det = J[0]*d[0*3 + 0] + J[3]*d[1*3 + 0] + J[6]*d[2*3 + 0];
}

//--- NEW Computation of Jacobian inverses ---

/// Compute Jacobian inverse K for interval embedded in R^1
inline void new_compute_jacobian_inverse_interval_1d(double K[UFC_TDIM_1*UFC_GDIM_1],
                                                     double det)
{
  K[0] = 1.0 / det;
}

/// Compute Jacobian (pseudo)inverse K for interval embedded in R^2
inline void new_compute_jacobian_inverse_interval_2d(double K[UFC_TDIM_1*UFC_GDIM_2],
                                                     double det2,
                                                     const double J[UFC_GDIM_2*UFC_TDIM_1])
{
  K[0] = J[0] / det2;
  K[1] = J[1] / det2;
}

/// Compute Jacobian (pseudo)inverse K for interval embedded in R^3
inline void new_compute_jacobian_inverse_interval_3d(double K[UFC_TDIM_1*UFC_GDIM_3],
                                                     double det2,
                                                     const double J[UFC_GDIM_3*UFC_TDIM_1])
{
  K[0] = J[0] / det2;
  K[1] = J[1] / det2;
  K[2] = J[2] / det2;
}

/// Compute Jacobian inverse K for triangle embedded in R^2
inline void new_compute_jacobian_inverse_triangle_2d(double K[UFC_TDIM_2*UFC_GDIM_2],
                                                     double det,
                                                     const double J[UFC_GDIM_2*UFC_TDIM_2])
{
  K[0] =  J[3] / det;
  K[1] = -J[1] / det;
  K[2] = -J[2] / det;
  K[3] =  J[0] / det;
}

/// Compute Jacobian (pseudo)inverse K for triangle embedded in R^3
inline void new_compute_jacobian_inverse_triangle_3d(double K[UFC_TDIM_2*UFC_GDIM_3],
                                                     double den,
                                                     const double c[3],
                                                     const double J[UFC_GDIM_3*UFC_TDIM_2])
{
  K[0] = (J[0]*c[1] - J[1]*c[2]) / den;
  K[1] = (J[2]*c[1] - J[3]*c[2]) / den;
  K[2] = (J[4]*c[1] - J[5]*c[2]) / den;
  K[3] = (J[1]*c[0] - J[0]*c[2]) / den;
  K[4] = (J[3]*c[0] - J[2]*c[2]) / den;
  K[5] = (J[5]*c[0] - J[4]*c[2]) / den;
}

/// Compute Jacobian inverse K for tetrahedron embedded in R^3
inline void new_compute_jacobian_inverse_tetrahedron_3d(double K[UFC_TDIM_3*UFC_GDIM_3],
                                                        double det,
                                                        const double d[9])
{
  K[0] = d[0*3 + 0] / det;
  K[1] = d[1*3 + 0] / det;
  K[2] = d[2*3 + 0] / det;
  K[3] = d[0*3 + 1] / det;
  K[4] = d[1*3 + 1] / det;
  K[5] = d[2*3 + 1] / det;
  K[6] = d[0*3 + 2] / det;
  K[7] = d[1*3 + 2] / det;
  K[8] = d[2*3 + 2] / det;
}

// --- Computation of edge, face, facet scaling factors

/// Compute edge scaling factors for triangle embedded in R^2
inline void compute_edge_scaling_factors_triangle_2d(double dx[2],
                                                     const double coordinate_dofs[6],
                                                     std::size_t facet)
{
  // Get vertices on edge
  const unsigned int v0 = triangle_facet_vertices[facet][0];
  const unsigned int v1 = triangle_facet_vertices[facet][1];

  // Compute scale factor (length of edge scaled by length of reference interval)
  dx[0] = coordinate_dofs[2*v1 + 0] - coordinate_dofs[2*v0 + 0];
  dx[1] = coordinate_dofs[2*v1 + 1] - coordinate_dofs[2*v0 + 1];
}

/// Compute facet scaling factor for triangle embedded in R^2
inline void compute_facet_scaling_factor_triangle_2d(double & det,
                                                     const double dx[2])
{
  det = std::sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
}

/// Compute edge scaling factors for triangle embedded in R^3
inline void compute_edge_scaling_factors_triangle_3d(double dx[3],
                                                     const double coordinate_dofs[9],
                                                     std::size_t facet)
{
  // Get vertices on edge
  const unsigned int v0 = triangle_facet_vertices[facet][0];
  const unsigned int v1 = triangle_facet_vertices[facet][1];

  // Compute scale factor (length of edge scaled by length of reference interval)
  dx[0] = coordinate_dofs[3*v1 + 0] - coordinate_dofs[3*v0 + 0];
  dx[1] = coordinate_dofs[3*v1 + 1] - coordinate_dofs[3*v0 + 1];
  dx[2] = coordinate_dofs[3*v1 + 2] - coordinate_dofs[3*v0 + 2];
}

/// Compute facet scaling factor for triangle embedded in R^3
inline void compute_facet_scaling_factor_triangle_3d(double & det,
                                                     const double dx[3])
{
  det = std::sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
}

/// Compute face scaling factors for tetrahedron embedded in R^3
inline void compute_face_scaling_factors_tetrahedron_3d(double a[3],
                                                        const double coordinate_dofs[12],
                                                        std::size_t facet)
{
  // Get vertices on face
  const unsigned int v0 = tetrahedron_facet_vertices[facet][0];
  const unsigned int v1 = tetrahedron_facet_vertices[facet][1];
  const unsigned int v2 = tetrahedron_facet_vertices[facet][2];

  // Compute scale factor (area of face scaled by area of reference triangle)
  a[0] = (coordinate_dofs[3*v0 + 1]*coordinate_dofs[3*v1 + 2]  +
          coordinate_dofs[3*v0 + 2]*coordinate_dofs[3*v2 + 1]  +
          coordinate_dofs[3*v1 + 1]*coordinate_dofs[3*v2 + 2]) -
         (coordinate_dofs[3*v2 + 1]*coordinate_dofs[3*v1 + 2]  +
          coordinate_dofs[3*v2 + 2]*coordinate_dofs[3*v0 + 1]  +
          coordinate_dofs[3*v1 + 1]*coordinate_dofs[3*v0 + 2]);

  a[1] = (coordinate_dofs[3*v0 + 2]*coordinate_dofs[3*v1 + 0]  +
          coordinate_dofs[3*v0 + 0]*coordinate_dofs[3*v2 + 2]  +
          coordinate_dofs[3*v1 + 2]*coordinate_dofs[3*v2 + 0]) -
         (coordinate_dofs[3*v2 + 2]*coordinate_dofs[3*v1 + 0]  +
          coordinate_dofs[3*v2 + 0]*coordinate_dofs[3*v0 + 2]  +
          coordinate_dofs[3*v1 + 2]*coordinate_dofs[3*v0 + 0]);

  a[2] = (coordinate_dofs[3*v0 + 0]*coordinate_dofs[3*v1 + 1]  +
          coordinate_dofs[3*v0 + 1]*coordinate_dofs[3*v2 + 0]  +
          coordinate_dofs[3*v1 + 0]*coordinate_dofs[3*v2 + 1]) -
         (coordinate_dofs[3*v2 + 0]*coordinate_dofs[3*v1 + 1]  +
          coordinate_dofs[3*v2 + 1]*coordinate_dofs[3*v0 + 0]  +
          coordinate_dofs[3*v1 + 0]*coordinate_dofs[3*v0 + 1]);
}

/// Compute facet scaling factor for tetrahedron embedded in R^3
inline void compute_facet_scaling_factor_tetrahedron_3d(double & det,
                                                        const double a[3])
{
  det = std::sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

///--- Compute facet normal directions ---

/// Compute facet direction for interval embedded in R^1
inline void compute_facet_normal_direction_interval_1d(bool & direction,
                                                       const double coordinate_dofs[2],
                                                       std::size_t facet)
{
  direction = facet == 0
    ? coordinate_dofs[0] > coordinate_dofs[1]
    : coordinate_dofs[1] > coordinate_dofs[0];
}

/// Compute facet direction for triangle embedded in R^2
inline void compute_facet_normal_direction_triangle_2d(bool & direction,
                                                       const double coordinate_dofs[6],
                                                       const double dx[2],
                                                       std::size_t facet)
{
  const unsigned int v0 = triangle_facet_vertices[facet][0];
  direction = dx[1]*(coordinate_dofs[2*facet    ] - coordinate_dofs[2*v0    ])
            - dx[0]*(coordinate_dofs[2*facet + 1] - coordinate_dofs[2*v0 + 1])
            < 0;
}

/// Compute facet direction for tetrahedron embedded in R^3
inline void compute_facet_normal_direction_tetrahedron_3d(bool & direction,
                                                          const double coordinate_dofs[9],
                                                          const double a[3],
                                                          std::size_t facet)
{
  const unsigned int v0 = tetrahedron_facet_vertices[facet][0];
  direction = a[0]*(coordinate_dofs[3*facet    ] - coordinate_dofs[3*v0    ])
            + a[1]*(coordinate_dofs[3*facet + 1] - coordinate_dofs[3*v0 + 1])
            + a[2]*(coordinate_dofs[3*facet + 2] - coordinate_dofs[3*v0 + 2])
            < 0;
}

///--- Compute facet normal vectors ---

/// Compute facet normal for interval embedded in R^1
inline void compute_facet_normal_interval_1d(double n[UFC_GDIM_1],
                                             bool direction)
{
  // Facet normals are 1.0 or -1.0:   (-1.0) <-- X------X --> (1.0)
  n[0] = direction ? 1.0 : -1.0;
}

/// Compute facet normal for interval embedded in R^2
inline void compute_facet_normal_interval_2d(double n[UFC_GDIM_2],
                                             const double coordinate_dofs[4],
                                             std::size_t facet)
{
  if (facet == 0)
  {
    n[0] = coordinate_dofs[0] - coordinate_dofs[2];
    n[1] = coordinate_dofs[1] - coordinate_dofs[3];
  }
  else
  {
    n[0] = coordinate_dofs[2] - coordinate_dofs[0];
    n[1] = coordinate_dofs[3] - coordinate_dofs[1];
  }
  const double n_length = std::sqrt(n[0]*n[0] + n[1]*n[1]);
  n[0] /= n_length;
  n[1] /= n_length;
}

/// Compute facet normal for interval embedded in R^3
inline void compute_facet_normal_interval_3d(double n[UFC_GDIM_3],
                                             const double coordinate_dofs[6],
                                             std::size_t facet)
{
  if (facet == 0)
  {
    n[0] = coordinate_dofs[0] - coordinate_dofs[3];
    n[1] = coordinate_dofs[1] - coordinate_dofs[4];
    n[1] = coordinate_dofs[2] - coordinate_dofs[5];
  }
  else
  {
    n[0] = coordinate_dofs[3] - coordinate_dofs[0];
    n[1] = coordinate_dofs[4] - coordinate_dofs[1];
    n[1] = coordinate_dofs[5] - coordinate_dofs[2];
  }
  const double n_length = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  n[0] /= n_length;
  n[1] /= n_length;
  n[2] /= n_length;
}

/// Compute facet normal for triangle embedded in R^2
inline void compute_facet_normal_triangle_2d(double n[UFC_GDIM_2],
                                             const double dx[2],
                                             const double det,
                                             bool direction)
{
  // Compute facet normals from the facet scale factor constants
  n[0] = direction ?  dx[1] / det : -dx[1] / det;
  n[1] = direction ? -dx[0] / det :  dx[0] / det;
}


/// Compute facet normal for triangle embedded in R^3
inline void compute_facet_normal_triangle_3d(double n[UFC_GDIM_3],
                                             const double coordinate_dofs[6],
                                             std::size_t facet)
{
  // Compute facet normal for triangles in 3D
  const unsigned int vertex0 = facet;

  // Get coordinates corresponding the vertex opposite this
  const unsigned int vertex1 = triangle_facet_vertices[facet][0];
  const unsigned int vertex2 = triangle_facet_vertices[facet][1];

  // Define vectors n = (p2 - p0) and t = normalized (p2 - p1)
  n[0] = coordinate_dofs[3*vertex2 + 0] - coordinate_dofs[3*vertex0 + 0];
  n[1] = coordinate_dofs[3*vertex2 + 1] - coordinate_dofs[3*vertex0 + 1];
  n[2] = coordinate_dofs[3*vertex2 + 2] - coordinate_dofs[3*vertex0 + 2];

  double t0 = coordinate_dofs[3*vertex2 + 0] - coordinate_dofs[3*vertex1 + 0];
  double t1 = coordinate_dofs[3*vertex2 + 1] - coordinate_dofs[3*vertex1 + 1];
  double t2 = coordinate_dofs[3*vertex2 + 2] - coordinate_dofs[3*vertex1 + 2];
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

/// Compute facet normal for tetrahedron embedded in R^3
inline void compute_facet_normal_tetrahedron_3d(double n[UFC_GDIM_3],
                                                const double a[3],
                                                const double det,
                                                bool direction)
{
  // Compute facet normals from the facet scale factor constants
  n[0] = direction ? a[0] / det : -a[0] / det;
  n[1] = direction ? a[1] / det : -a[1] / det;
  n[2] = direction ? a[2] / det : -a[2] / det;
}

///--- Compute circumradius ---

/// Compute circumradius for interval embedded in R^1
inline void compute_circumradius_interval_1d(double & circumradius,
                                             double volume)
{
  // Compute circumradius; in 1D it is equal to half the cell length
  circumradius = volume / 2.0;
}


/// Compute circumradius for interval embedded in R^2
inline void compute_circumradius_interval_2d(double & circumradius,
                                             double volume)
{
  // Compute circumradius of interval in 2D (1/2 volume)
  circumradius = volume / 2.0;
}


/// Compute circumradius for interval embedded in R^3
inline void compute_circumradius_interval_3d(double & circumradius,
                                             double volume)
{
  // Compute circumradius of interval in 3D (1/2 volume)
  circumradius = volume / 2.0;
}

/// Compute circumradius for triangle embedded in R^2
inline void compute_circumradius_triangle_2d(double & circumradius,
                                             const double coordinate_dofs[6],
                                             const double J[UFC_GDIM_2*UFC_TDIM_2],
                                             double volume)
{
  // Compute circumradius of triangle in 2D
  const double v1v2  = std::sqrt(  (coordinate_dofs[4] - coordinate_dofs[2])*(coordinate_dofs[4] - coordinate_dofs[2])
                                 + (coordinate_dofs[5] - coordinate_dofs[3])*(coordinate_dofs[5] - coordinate_dofs[3]) );
  const double v0v2  = std::sqrt(J[3]*J[3] + J[1]*J[1]);
  const double v0v1  = std::sqrt(J[0]*J[0] + J[2]*J[2]);

  circumradius = 0.25*(v1v2*v0v2*v0v1) / volume;
}

/// Compute circumradius for triangle embedded in R^3
inline void compute_circumradius_triangle_3d(double & circumradius,
                                             const double coordinate_dofs[9],
                                             const double J[UFC_GDIM_3*UFC_TDIM_2],
                                             double volume)
{
  // Compute circumradius of triangle in 3D
  const double v1v2  = std::sqrt(   (coordinate_dofs[6] - coordinate_dofs[3])*(coordinate_dofs[6] - coordinate_dofs[3])
                                  + (coordinate_dofs[7] - coordinate_dofs[4])*(coordinate_dofs[7] - coordinate_dofs[4])
                                  + (coordinate_dofs[8] - coordinate_dofs[5])*(coordinate_dofs[8] - coordinate_dofs[5]));
  const double v0v2 = std::sqrt( J[3]*J[3] + J[1]*J[1] + J[5]*J[5]);
  const double v0v1 = std::sqrt( J[0]*J[0] + J[2]*J[2] + J[4]*J[4]);

  circumradius = 0.25*(v1v2*v0v2*v0v1) / volume;
}

/// Compute circumradius for tetrahedron embedded in R^3
inline void compute_circumradius_tetrahedron_3d(double & circumradius,
                                                const double coordinate_dofs[12],
                                                const double J[UFC_GDIM_3*UFC_TDIM_3],
                                                double volume)
{
  // Compute circumradius
  const double v1v2  = std::sqrt(   (coordinate_dofs[6] - coordinate_dofs[3])*(coordinate_dofs[6] - coordinate_dofs[3])
                                  + (coordinate_dofs[7] - coordinate_dofs[4])*(coordinate_dofs[7] - coordinate_dofs[4])
                                  + (coordinate_dofs[8] - coordinate_dofs[5])*(coordinate_dofs[8] - coordinate_dofs[5]) );
  const double v0v2  = std::sqrt(J[1]*J[1] + J[4]*J[4] + J[7]*J[7]);
  const double v0v1  = std::sqrt(J[0]*J[0] + J[3]*J[3] + J[6]*J[6]);
  const double v0v3  = std::sqrt(J[2]*J[2] + J[5]*J[5] + J[8]*J[8]);
  const double v1v3  = std::sqrt(   (coordinate_dofs[ 9] - coordinate_dofs[3])*(coordinate_dofs[ 9] - coordinate_dofs[3])
                                  + (coordinate_dofs[10] - coordinate_dofs[4])*(coordinate_dofs[10] - coordinate_dofs[4])
                                  + (coordinate_dofs[11] - coordinate_dofs[5])*(coordinate_dofs[11] - coordinate_dofs[5]) );
  const double v2v3  = std::sqrt(   (coordinate_dofs[ 9] - coordinate_dofs[6])*(coordinate_dofs[ 9] - coordinate_dofs[6])
                                  + (coordinate_dofs[10] - coordinate_dofs[7])*(coordinate_dofs[10] - coordinate_dofs[7])
                                  + (coordinate_dofs[11] - coordinate_dofs[8])*(coordinate_dofs[11] - coordinate_dofs[8]) );
  const  double la   = v1v2*v0v3;
  const  double lb   = v0v2*v1v3;
  const  double lc   = v0v1*v2v3;
  const  double s    = 0.5*(la+lb+lc);
  const  double area = std::sqrt(s*(s-la)*(s-lb)*(s-lc));

  circumradius = area / (6.0*volume);
}

///--- Compute max facet edge lengths ---

/// Compute min edge length in facet of tetrahedron embedded in R^3
inline void compute_min_facet_edge_length_tetrahedron_3d(double & min_edge_length,
                                                         unsigned int facet,
                                                         const double coordinate_dofs[3*4])
{
  // TODO: Extract compute_facet_edge_lengths_tetrahedron_3d(), reuse between min/max functions
  double edge_lengths_sqr[3];
  for (unsigned int edge = 0; edge < 3; ++edge)
  {
    const unsigned int vertex0 = tetrahedron_facet_edge_vertices[facet][edge][0];
    const unsigned int vertex1 = tetrahedron_facet_edge_vertices[facet][edge][1];
    edge_lengths_sqr[edge] = (coordinate_dofs[3*vertex1 + 0] - coordinate_dofs[3*vertex0 + 0])*(coordinate_dofs[3*vertex1 + 0] - coordinate_dofs[3*vertex0 + 0])
                           + (coordinate_dofs[3*vertex1 + 1] - coordinate_dofs[3*vertex0 + 1])*(coordinate_dofs[3*vertex1 + 1] - coordinate_dofs[3*vertex0 + 1])
                           + (coordinate_dofs[3*vertex1 + 2] - coordinate_dofs[3*vertex0 + 2])*(coordinate_dofs[3*vertex1 + 2] - coordinate_dofs[3*vertex0 + 2]);
  }
  min_edge_length = std::sqrt(std::min(std::min(edge_lengths_sqr[1], edge_lengths_sqr[1]), edge_lengths_sqr[2]));
}

///--- Compute max facet edge lengths ---

/// Compute max edge length in facet of tetrahedron embedded in R^3
inline void compute_max_facet_edge_length_tetrahedron_3d(double & max_edge_length,
                                                         unsigned int facet,
                                                         const double coordinate_dofs[12])
{
  // TODO: Extract compute_facet_edge_lengths_tetrahedron_3d(), reuse between min/max functions
  double edge_lengths_sqr[3];
  for (unsigned int edge = 0; edge < 3; ++edge)
  {
    const unsigned int vertex0 = tetrahedron_facet_edge_vertices[facet][edge][0];
    const unsigned int vertex1 = tetrahedron_facet_edge_vertices[facet][edge][1];
    edge_lengths_sqr[edge] = (coordinate_dofs[3*vertex1 + 0] - coordinate_dofs[3*vertex0 + 0])*(coordinate_dofs[3*vertex1 + 0] - coordinate_dofs[3*vertex0 + 0])
                           + (coordinate_dofs[3*vertex1 + 1] - coordinate_dofs[3*vertex0 + 1])*(coordinate_dofs[3*vertex1 + 1] - coordinate_dofs[3*vertex0 + 1])
                           + (coordinate_dofs[3*vertex1 + 2] - coordinate_dofs[3*vertex0 + 2])*(coordinate_dofs[3*vertex1 + 2] - coordinate_dofs[3*vertex0 + 2]);
  }
  max_edge_length = std::sqrt(std::max(std::max(edge_lengths_sqr[0], edge_lengths_sqr[1]), edge_lengths_sqr[2]));
}


//} // TODO: Wrap all in namespace ufc


#endif
