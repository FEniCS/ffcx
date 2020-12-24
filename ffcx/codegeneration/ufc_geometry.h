// This file provides utility functions for computing geometric quantities.
// This code is released into the public domain.
//
// The FEniCS Project (http://www.fenicsproject.org/) 2013-2018.

#pragma once

#include <math.h>
#include <stdbool.h>

/// A note regarding data structures. All matrices are represented as
/// row-major flattened raw C arrays.

/// TODO: Remove this file and get data from Basix instead

// TODO: Should signatures of compute_<foo>_<cell>_<n>d match for each foo?
//       On one hand the snippets use different quantities, on the other
//       some consistency is nice to simplify the code generation.
//       Currently only the arguments that are actually used are included.

/// --- Some fixed numbers by name for readability in this file ---
// TODO: Use these numbers where relevant below to make this file more self
// documenting

// (namespaced using UFC_ in the names because they collide with variables in
// other libraries)

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

/// --- Local reference cell coordinates by basix conventions ---

static const double interval_vertices[UFC_NUM_VERTICES_IN_INTERVAL][UFC_TDIM_1]
    = {{0.0}, {1.0}};

static const double triangle_vertices[UFC_NUM_VERTICES_IN_TRIANGLE][UFC_TDIM_2]
    = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}};

static const double tetrahedron_vertices[UFC_NUM_VERTICES_IN_TETRAHEDRON]
                                        [UFC_TDIM_3]
    = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

static const double quadrilateral_vertices[UFC_NUM_VERTICES_IN_QUADRILATERAL]
                                          [UFC_TDIM_2]
    = {
        {0.0, 0.0},
        {1.0, 0.0},
        {0.0, 1.0},
        {1.0, 1.0},
};

static const double hexahedron_vertices[UFC_NUM_VERTICES_IN_HEXAHEDRON]
                                       [UFC_TDIM_3]
    = {
        {0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 1.0, 0.0}, {0.0, 1.0, 1.0},
        {1.0, 0.0, 0.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 0.0}, {1.0, 1.0, 1.0},
};

/// --- Local reference cell midpoint by basix conventions ---

static const double interval_midpoint[UFC_TDIM_1] = {0.5};

static const double triangle_midpoint[UFC_TDIM_2] = {1.0 / 3.0, 1.0 / 3.0};

static const double tetrahedron_midpoint[UFC_TDIM_3] = {0.25, 0.25, 0.25};

static const double quadrilateral_midpoint[UFC_TDIM_2] = {0.5, 0.5};

static const double hexahedron_midpoint[UFC_TDIM_3] = {0.5, 0.5, 0.5};

/// --- Local reference cell facet midpoints by basix conventions ---

static const double interval_facet_midpoint[UFC_NUM_FACETS_IN_INTERVAL]
                                           [UFC_TDIM_1]
    = {{0.0}, {1.0}};

static const double triangle_facet_midpoint[UFC_NUM_FACETS_IN_TRIANGLE]
                                           [UFC_TDIM_2]
    = {{0.5, 0.5}, {0.0, 0.5}, {0.5, 0.0}};

static const double tetrahedron_facet_midpoint[UFC_NUM_FACETS_IN_TETRAHEDRON]
                                              [UFC_TDIM_3]
    = {
        {0.5, 0.5, 0.5},
        {0.0, 1.0 / 3.0, 1.0 / 3.0},
        {1.0 / 3.0, 0.0, 1.0 / 3.0},
        {1.0 / 3.0, 1.0 / 3.0, 0.0},
};

static const double
    quadrilateral_facet_midpoint[UFC_NUM_FACETS_IN_QUADRILATERAL][UFC_TDIM_2]
    = {
        {0.5, 0.0},
        {0.0, 0.5},
        {1.0, 0.5},
        {0.5, 1.0},
};

static const double hexahedron_facet_midpoint[UFC_NUM_FACETS_IN_HEXAHEDRON]
                                             [UFC_TDIM_3]
    = {
        {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5},
        {1.0, 0.5, 0.5}, {0.5, 1.0, 0.5}, {0.5, 0.5, 1.0},
};

/// --- Local reference cell facet orientations by basix conventions ---

static const double interval_facet_orientations[UFC_NUM_FACETS_IN_INTERVAL] = {
    -1.0,
    +1.0,
};

static const double triangle_facet_orientations[UFC_NUM_FACETS_IN_TRIANGLE]
    = {+1.0, -1.0, +1.0};

static const double
    tetrahedron_facet_orientations[UFC_NUM_FACETS_IN_TETRAHEDRON]
    = {+1.0, -1.0, +1.0, -1.0};

static const double
quadrilateral_facet_orientations[UFC_NUM_FACETS_IN_QUADRILATERAL] = {
  -1.0,
  +1.0,
  -1.0,
  +1.0,
};

static const double hexahedron_facet_orientations[UFC_NUM_FACETS_IN_HEXAHEDRON]
= {
  -1.0,
  +1.0,
  -1.0,
  +1.0,
  -1.0,
  +1.0,
};

/// --- Local reference cell entity relations by basix conventions ---

static const unsigned int triangle_edge_vertices[UFC_NUM_EDGES_IN_TRIANGLE][2]
    = {{1, 2}, {0, 2}, {0, 1}};

static const unsigned int
    tetrahedron_edge_vertices[UFC_NUM_EDGES_IN_TETRAHEDRON][2]
    = {{2, 3}, {1, 3}, {1, 2}, {0, 3}, {0, 2}, {0, 1}};

static const unsigned int
    quadrilateral_edge_vertices[UFC_NUM_EDGES_IN_QUADRILATERAL][2]
    = {
        {0, 1},
        {0, 2},
        {1, 3},
        {2, 3},
};

static const unsigned int hexahedron_edge_vertices[UFC_NUM_EDGES_IN_HEXAHEDRON]
                                                  [2]
    = {
        {0, 1}, {0, 2}, {0, 4}, {1, 3}, {1, 5}, {2, 3},
        {2, 6}, {3, 7}, {4, 5}, {4, 7}, {5, 7}, {6, 7},
};

/// --- Local reference cell entity relations by basix conventions ---

static const unsigned int interval_facet_vertices[UFC_NUM_FACETS_IN_INTERVAL][1]
    = {{0}, {1}};

static const unsigned int triangle_facet_vertices[UFC_NUM_FACETS_IN_TRIANGLE]
                                                 [UFC_NUM_VERTICES_IN_INTERVAL]
    = {{1, 2}, {0, 2}, {0, 1}};

static const unsigned int
    tetrahedron_facet_vertices[UFC_NUM_FACETS_IN_TETRAHEDRON]
                              [UFC_NUM_VERTICES_IN_TRIANGLE]
    = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};

static const unsigned int
    tetrahedron_facet_edge_vertices[UFC_NUM_FACETS_IN_TETRAHEDRON]
                                   [UFC_NUM_FACETS_IN_TRIANGLE]
                                   [UFC_NUM_VERTICES_IN_INTERVAL]
    = {{{2, 3}, {1, 3}, {1, 2}},
       {{2, 3}, {0, 3}, {0, 2}},
       {{1, 3}, {0, 3}, {0, 1}},
       {{1, 2}, {0, 2}, {0, 1}}};

static const unsigned int
    quadrilateral_facet_vertices[UFC_NUM_FACETS_IN_QUADRILATERAL]
                                [UFC_NUM_VERTICES_IN_INTERVAL]
    = {
        {0, 1},
        {0, 2},
        {1, 3},
        {2, 3},
};

static const unsigned int
    hexahedron_facet_vertices[UFC_NUM_FACETS_IN_HEXAHEDRON]
                             [UFC_NUM_VERTICES_IN_QUADRILATERAL]
    = {
        {0, 1, 2, 3}, {0, 1, 4, 5}, {0, 2, 4, 6},
        {1, 3, 5, 7}, {2, 3, 6, 7}, {4, 5, 6, 7},
};

static const unsigned int
    hexahedron_facet_edge_vertices[UFC_NUM_FACETS_IN_HEXAHEDRON]
                                  [UFC_NUM_FACETS_IN_QUADRILATERAL]
                                  [UFC_NUM_VERTICES_IN_INTERVAL]
    = {
        {{0, 1}, {0, 2}, {1, 3}, {2, 3}}, {{0, 1}, {0, 4}, {1, 5}, {4, 5}},
        {{0, 2}, {0, 4}, {2, 6}, {4, 6}}, {{1, 3}, {1, 5}, {3, 7}, {5, 7}},
        {{2, 3}, {2, 6}, {3, 6}, {3, 7}}, {{4, 5}, {4, 6}, {5, 7}, {6, 7}},
};

/// --- Reference cell edge vectors by basix conventions (edge vertex 1 - edge
/// vertex 0 for each edge in cell) ---

static const double triangle_reference_edge_vectors[UFC_NUM_EDGES_IN_TRIANGLE]
                                                   [UFC_TDIM_2]
    = {
        {-1.0, 1.0},
        {0.0, 1.0},
        {1.0, 0.0},
};

static const double
    tetrahedron_reference_edge_vectors[UFC_NUM_EDGES_IN_TETRAHEDRON][UFC_TDIM_3]
    = {
        {0.0, -1.0, 1.0}, {-1.0, 0.0, 1.0}, {-1.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},  {0.0, 1.0, 0.0},  {1.0, 0.0, 0.0},
};

// Edge vectors for each triangle facet of a tetrahedron
static const double tetrahedron_facet_reference_edge_vectors
    [UFC_NUM_FACETS_IN_TETRAHEDRON][UFC_NUM_EDGES_IN_TRIANGLE][UFC_TDIM_3]
    = {
        {
            // facet 0
            {0.0, -1.0, 1.0},
            {-1.0, 0.0, 1.0},
            {-1.0, 1.0, 0.0},
        },
        {
            // facet 1
            {0.0, -1.0, 1.0},
            {0.0, 0.0, 1.0},
            {0.0, 1.0, 0.0},
        },
        {
            // facet 2
            {-1.0, 0.0, 1.0},
            {0.0, 0.0, 1.0},
            {1.0, 0.0, 0.0},
        },
        {
            // facet 3
            {-1.0, 1.0, 0.0},
            {0.0, 1.0, 0.0},
            {1.0, 0.0, 0.0},
        },
};

static const double
    quadrilateral_reference_edge_vectors[UFC_NUM_EDGES_IN_QUADRILATERAL]
                                        [UFC_TDIM_2]
    = {
        {1.0, 0.0},
        {0.0, 1.0},
        {0.0, 1.0},
        {1.0, 0.0},
};

static const double
    hexahedron_reference_edge_vectors[UFC_NUM_EDGES_IN_HEXAHEDRON][UFC_TDIM_3]
    = {
        {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 0.0, 1.0},
        {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 0.0, 0.0},
};

// Edge vectors for each quadrilateral facet of a hexahedron
static const double hexahedron_facet_reference_edge_vectors
    [UFC_NUM_FACETS_IN_HEXAHEDRON][UFC_NUM_EDGES_IN_QUADRILATERAL][UFC_TDIM_3]
    = {
        {
            // facet 0
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 1.0, 0.0},
            {1.0, 0.0, 0.0},
        },
        {
            // facet 1
            {1.0, 0.0, 0.0},
            {0.0, 0.0, 1.0},
            {0.0, 0.0, 1.0},
            {1.0, 0.0, 0.0},
        },
        {
            // facet 2
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0},
            {0.0, 0.0, 1.0},
            {0.0, 1.0, 0.0},
        },
        {
            // facet 3
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0},
            {0.0, 0.0, 1.0},
            {0.0, 1.0, 0.0},
        },
        {
            // facet 4
            {1.0, 0.0, 0.0},
            {0.0, 0.0, 1.0},
            {0.0, 0.0, 1.0},
            {1.0, 0.0, 0.0},
        },
        {
            // facet 5
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 1.0, 0.0},
            {1.0, 0.0, 0.0},
        },
};

/// --- Reference cell facet normals by basix conventions (outwards pointing on
/// reference cell) ---

static const double interval_reference_facet_normals[UFC_NUM_FACETS_IN_INTERVAL]
                                                    [UFC_TDIM_1]
    = {
        {-1.0},
        {+1.0},
};

static const double triangle_reference_facet_normals[UFC_NUM_FACETS_IN_TRIANGLE]
                                                    [UFC_TDIM_2]
    = {
        {0.7071067811865476, 0.7071067811865476},
        {-1.0, 0.0},
        {0.0, -1.0},
};

static const double
    tetrahedron_reference_facet_normals[UFC_NUM_FACETS_IN_TETRAHEDRON]
                                       [UFC_TDIM_3]
    = {
        {0.5773502691896258, 0.5773502691896258, 0.5773502691896258},
        {-1.0, 0.0, 0.0},
        {0.0, -1.0, 0.0},
        {0.0, 0.0, -1.0},
};

static const double
    quadrilateral_reference_facet_normals[UFC_NUM_FACETS_IN_QUADRILATERAL]
                                         [UFC_TDIM_2]
    = {
        {0.0, 1.0},
        {-1.0, 0.0},
        {-1.0, 0.0},
        {0.0, 1.0},
};

static const double
    hexahedron_reference_facet_normals[UFC_NUM_FACETS_IN_HEXAHEDRON][UFC_TDIM_3]
    = {
        {0.0, 0.0, -1.0}, {0.0, -1.0, 0.0},  {-1.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},  {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0},
};

/// --- Reference cell volumes by basix conventions ---

static const double interval_reference_cell_volume = 1.0;
static const double triangle_reference_cell_volume = 0.5;
static const double tetrahedron_reference_cell_volume = 1.0 / 6.0;
static const double quadrilateral_reference_cell_volume = 1.0;
static const double hexahedron_reference_cell_volume = 1.0;

static const double interval_reference_facet_volume = 1.0;
static const double triangle_reference_facet_volume = 1.0;
static const double tetrahedron_reference_facet_volume = 0.5;
static const double quadrilateral_reference_facet_volume = 1.0;
static const double hexahedron_reference_facet_volume = 1.0;

/// --- Jacobians of reference facet cell to reference cell coordinate mappings
/// by basix conventions ---

static const double
    triangle_reference_facet_jacobian[UFC_NUM_FACETS_IN_TRIANGLE][UFC_TDIM_2]
                                     [UFC_TDIM_2 - 1]
    = {
        {{-1.0}, {1.0}},
        {{0.0}, {1.0}},
        {{1.0}, {0.0}},
};

static const double
    tetrahedron_reference_facet_jacobian[UFC_NUM_FACETS_IN_TETRAHEDRON]
                                        [UFC_TDIM_3][UFC_TDIM_3 - 1]
    = {
        {{-1.0, -1.0}, {1.0, 0.0}, {0.0, 1.0}},
        {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}},
        {{1.0, 0.0}, {0.0, 0.0}, {0.0, 1.0}},
        {{1.0, 0.0}, {0.0, 1.0}, {0.0, 0.0}},
};

static const double
    quadrilateral_reference_facet_jacobian[UFC_NUM_FACETS_IN_QUADRILATERAL]
                                          [UFC_TDIM_2][UFC_TDIM_2 - 1]
    = {
        {{1.0}, {0.0}},
        {{0.0}, {1.0}},
        {{0.0}, {1.0}},
        {{1.0}, {0.0}},
};

static const double
    hexahedron_reference_facet_jacobian[UFC_NUM_FACETS_IN_HEXAHEDRON]
                                       [UFC_TDIM_3][UFC_TDIM_3 - 1]
    = {
        {{1.0, 0.0}, {0.0, 1.0}, {0.0, 0.0}},
        {{1.0, 0.0}, {0.0, 0.0}, {0.0, 1.0}},
        {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}},
        {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}},
        {{1.0, 0.0}, {0.0, 0.0}, {0.0, 1.0}},
        {{1.0, 0.0}, {0.0, 1.0}, {0.0, 0.0}},
};
