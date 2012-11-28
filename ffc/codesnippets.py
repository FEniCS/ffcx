"Code snippets for code generation."

# Copyright (C) 2007 Anders Logg
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Kristian B. Oelgaard 2010-2011
# Modified by Marie Rognes 2007-2010
# Modified by Peter Brune 2009
#
# First added:  2007-02-28
# Last changed: 2011-11-22

# Code snippets

__all__ = ["comment_ufc", "comment_dolfin", "header_h", "header_c", "footer",
           "cell_coordinates", "jacobian", "inverse_jacobian",
           "evaluate_f",
           "facet_determinant", "map_onto_physical",
           "fiat_coordinate_map", "transform_snippet",
           "scale_factor", "combinations_snippet",
           "normal_direction",
           "facet_normal", "ip_coordinates", "cell_volume", "circumradius",
           "facet_area"]

comment_ufc = """\
// This code conforms with the UFC specification version %(ufc_version)s
// and was automatically generated by FFC version %(ffc_version)s.
"""

comment_dolfin = """\
// This code conforms with the UFC specification version %(ufc_version)s
// and was automatically generated by FFC version %(ffc_version)s.
//
// This code was generated with the option '-l dolfin' and
// contains DOLFIN-specific wrappers that depend on DOLFIN.
"""

header_h = """\
#ifndef __%(prefix_upper)s_H
#define __%(prefix_upper)s_H

#include <cmath>
#include <stdexcept>
#include <fstream>
#include <ufc.h>
"""

header_c = """\
#include "%(prefix)s.h"
"""

footer = """\
#endif
"""

cell_coordinates = "const double * const * x = c.coordinates;\n"

# Code snippets for computing Jacobian
_jacobian_1D = """\
// Extract vertex coordinates
const double * const * x%(restriction)s = c%(restriction)s.coordinates;

// Compute Jacobian of affine map from reference cell
const double J%(restriction)s_00 = x%(restriction)s[1][0] - x%(restriction)s[0][0];"""

_jacobian_2D = """\
// Extract vertex coordinates
const double * const * x%(restriction)s = c%(restriction)s.coordinates;

// Compute Jacobian of affine map from reference cell
const double J%(restriction)s_00 = x%(restriction)s[1][0] - x%(restriction)s[0][0];
const double J%(restriction)s_01 = x%(restriction)s[2][0] - x%(restriction)s[0][0];
const double J%(restriction)s_10 = x%(restriction)s[1][1] - x%(restriction)s[0][1];
const double J%(restriction)s_11 = x%(restriction)s[2][1] - x%(restriction)s[0][1];"""

_jacobian_2D_1D = """\
// Geometric dimension 2, topological dimension 1

// Extract vertex coordinates
const double * const * x%(restriction)s = c%(restriction)s.coordinates;

// Compute Jacobian of affine map from reference cell
const double J%(restriction)s_00 = x%(restriction)s[1][0] - x%(restriction)s[0][0];
const double J%(restriction)s_10 = x%(restriction)s[1][1] - x%(restriction)s[0][1];
"""

_jacobian_3D = """\
// Extract vertex coordinates
const double * const * x%(restriction)s = c%(restriction)s.coordinates;

// Compute Jacobian of affine map from reference cell
const double J%(restriction)s_00 = x%(restriction)s[1][0] - x%(restriction)s[0][0];
const double J%(restriction)s_01 = x%(restriction)s[2][0] - x%(restriction)s[0][0];
const double J%(restriction)s_02 = x%(restriction)s[3][0] - x%(restriction)s[0][0];
const double J%(restriction)s_10 = x%(restriction)s[1][1] - x%(restriction)s[0][1];
const double J%(restriction)s_11 = x%(restriction)s[2][1] - x%(restriction)s[0][1];
const double J%(restriction)s_12 = x%(restriction)s[3][1] - x%(restriction)s[0][1];
const double J%(restriction)s_20 = x%(restriction)s[1][2] - x%(restriction)s[0][2];
const double J%(restriction)s_21 = x%(restriction)s[2][2] - x%(restriction)s[0][2];
const double J%(restriction)s_22 = x%(restriction)s[3][2] - x%(restriction)s[0][2];"""

_jacobian_3D_2D = """\
// Geometric dimension 3, topological dimension 2

// Extract vertex coordinates
const double * const * x%(restriction)s = c%(restriction)s.coordinates;

// Compute Jacobian of affine map from reference cell
const double J%(restriction)s_00 = x%(restriction)s[1][0] - x%(restriction)s[0][0];
const double J%(restriction)s_01 = x%(restriction)s[2][0] - x%(restriction)s[0][0];
const double J%(restriction)s_10 = x%(restriction)s[1][1] - x%(restriction)s[0][1];
const double J%(restriction)s_11 = x%(restriction)s[2][1] - x%(restriction)s[0][1];
const double J%(restriction)s_20 = x%(restriction)s[1][2] - x%(restriction)s[0][2];
const double J%(restriction)s_21 = x%(restriction)s[2][2] - x%(restriction)s[0][2];"""

_jacobian_3D_1D = """\
// Geometric dimension 3, topological dimension 1

// Extract vertex coordinates
const double * const * x%(restriction)s = c%(restriction)s.coordinates;

// Compute Jacobian of affine map from reference cell
const double J%(restriction)s_00 = x%(restriction)s[1][0] - x%(restriction)s[0][0];
const double J%(restriction)s_10 = x%(restriction)s[1][1] - x%(restriction)s[0][1];
const double J%(restriction)s_20 = x%(restriction)s[1][2] - x%(restriction)s[0][2];"""

# Code snippets for computing the inverse Jacobian. Assumes that
# Jacobian is already initialized
_inverse_jacobian_1D = """\

// Compute determinant of Jacobian
const double detJ%(restriction)s =  J%(restriction)s_00;

// Compute inverse of Jacobian
const double K%(restriction)s_00 =  1.0 / detJ%(restriction)s;"""

_inverse_jacobian_2D = """\

// Compute determinant of Jacobian
double detJ%(restriction)s = J%(restriction)s_00*J%(restriction)s_11 - J%(restriction)s_01*J%(restriction)s_10;

// Compute inverse of Jacobian
const double K%(restriction)s_00 =  J%(restriction)s_11 / detJ%(restriction)s;
const double K%(restriction)s_01 = -J%(restriction)s_01 / detJ%(restriction)s;
const double K%(restriction)s_10 = -J%(restriction)s_10 / detJ%(restriction)s;
const double K%(restriction)s_11 =  J%(restriction)s_00 / detJ%(restriction)s;"""

_inverse_jacobian_2D_1D = """\

// Compute determinant of Jacobian
double detJ%(restriction)s = std::sqrt(J%(restriction)s_00*J%(restriction)s_00 + J%(restriction)s_10*J%(restriction)s_10);

// Compute inverse of Jacobian
const double K%(restriction)s_00 =  J%(restriction)s_00 / std::pow(detJ%(restriction)s, 2);
const double K%(restriction)s_01 =  J%(restriction)s_10 / std::pow(detJ%(restriction)s, 2);"""

_inverse_jacobian_3D = """\

// Compute sub determinants
const double d%(restriction)s_00 = J%(restriction)s_11*J%(restriction)s_22 - J%(restriction)s_12*J%(restriction)s_21;
const double d%(restriction)s_01 = J%(restriction)s_12*J%(restriction)s_20 - J%(restriction)s_10*J%(restriction)s_22;
const double d%(restriction)s_02 = J%(restriction)s_10*J%(restriction)s_21 - J%(restriction)s_11*J%(restriction)s_20;
const double d%(restriction)s_10 = J%(restriction)s_02*J%(restriction)s_21 - J%(restriction)s_01*J%(restriction)s_22;
const double d%(restriction)s_11 = J%(restriction)s_00*J%(restriction)s_22 - J%(restriction)s_02*J%(restriction)s_20;
const double d%(restriction)s_12 = J%(restriction)s_01*J%(restriction)s_20 - J%(restriction)s_00*J%(restriction)s_21;
const double d%(restriction)s_20 = J%(restriction)s_01*J%(restriction)s_12 - J%(restriction)s_02*J%(restriction)s_11;
const double d%(restriction)s_21 = J%(restriction)s_02*J%(restriction)s_10 - J%(restriction)s_00*J%(restriction)s_12;
const double d%(restriction)s_22 = J%(restriction)s_00*J%(restriction)s_11 - J%(restriction)s_01*J%(restriction)s_10;

// Compute determinant of Jacobian
double detJ%(restriction)s = J%(restriction)s_00*d%(restriction)s_00 + J%(restriction)s_10*d%(restriction)s_10 + J%(restriction)s_20*d%(restriction)s_20;

// Compute inverse of Jacobian
const double K%(restriction)s_00 = d%(restriction)s_00 / detJ%(restriction)s;
const double K%(restriction)s_01 = d%(restriction)s_10 / detJ%(restriction)s;
const double K%(restriction)s_02 = d%(restriction)s_20 / detJ%(restriction)s;
const double K%(restriction)s_10 = d%(restriction)s_01 / detJ%(restriction)s;
const double K%(restriction)s_11 = d%(restriction)s_11 / detJ%(restriction)s;
const double K%(restriction)s_12 = d%(restriction)s_21 / detJ%(restriction)s;
const double K%(restriction)s_20 = d%(restriction)s_02 / detJ%(restriction)s;
const double K%(restriction)s_21 = d%(restriction)s_12 / detJ%(restriction)s;
const double K%(restriction)s_22 = d%(restriction)s_22 / detJ%(restriction)s;"""

evaluate_f = "f.evaluate(vals, y, c);"

scale_factor = """\
// Set scale factor
const double det = std::abs(detJ);"""

_facet_determinant_1D = """\
// Facet determinant 1D (vertex)
const double det = 1.0;"""

_facet_determinant_2D = """\
// Get vertices on edge
static unsigned int edge_vertices[3][2] = {{1, 2}, {0, 2}, {0, 1}};
const unsigned int v0 = edge_vertices[facet%(restriction)s][0];
const unsigned int v1 = edge_vertices[facet%(restriction)s][1];

// Compute scale factor (length of edge scaled by length of reference interval)
const double dx0 = x%(restriction)s[v1][0] - x%(restriction)s[v0][0];
const double dx1 = x%(restriction)s[v1][1] - x%(restriction)s[v0][1];
const double det = std::sqrt(dx0*dx0 + dx1*dx1);"""

_facet_determinant_3D = """\
// Get vertices on face
static unsigned int face_vertices[4][3] = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};
const unsigned int v0 = face_vertices[facet%(restriction)s][0];
const unsigned int v1 = face_vertices[facet%(restriction)s][1];
const unsigned int v2 = face_vertices[facet%(restriction)s][2];

// Compute scale factor (area of face scaled by area of reference triangle)
const double a0 = (x%(restriction)s[v0][1]*x%(restriction)s[v1][2] + x%(restriction)s[v0][2]*x%(restriction)s[v2][1] + x%(restriction)s[v1][1]*x%(restriction)s[v2][2]) - (x%(restriction)s[v2][1]*x%(restriction)s[v1][2] + x%(restriction)s[v2][2]*x%(restriction)s[v0][1] + x%(restriction)s[v1][1]*x%(restriction)s[v0][2]);

const double a1 = (x%(restriction)s[v0][2]*x%(restriction)s[v1][0] + x%(restriction)s[v0][0]*x%(restriction)s[v2][2] + x%(restriction)s[v1][2]*x%(restriction)s[v2][0]) - (x%(restriction)s[v2][2]*x%(restriction)s[v1][0] + x%(restriction)s[v2][0]*x%(restriction)s[v0][2] + x%(restriction)s[v1][2]*x%(restriction)s[v0][0]);

const double a2 = (x%(restriction)s[v0][0]*x%(restriction)s[v1][1] + x%(restriction)s[v0][1]*x%(restriction)s[v2][0] + x%(restriction)s[v1][0]*x%(restriction)s[v2][1]) - (x%(restriction)s[v2][0]*x%(restriction)s[v1][1] + x%(restriction)s[v2][1]*x%(restriction)s[v0][0] + x%(restriction)s[v1][0]*x%(restriction)s[v0][1]);

const double det = std::sqrt(a0*a0 + a1*a1 + a2*a2);"""

_normal_direction_1D = """\
const bool direction = facet%(restriction)s == 0 ? x%(restriction)s[0][0] > x%(restriction)s[1][0] : x%(restriction)s[1][0] > x%(restriction)s[0][0];"""

_normal_direction_2D = """\
const bool direction = dx1*(x%(restriction)s[%(facet)s][0] - x%(restriction)s[v0][0]) - dx0*(x%(restriction)s[%(facet)s][1] - x%(restriction)s[v0][1]) < 0;"""

_normal_direction_3D = """\
const bool direction = a0*(x%(restriction)s[%(facet)s][0] - x%(restriction)s[v0][0]) + a1*(x%(restriction)s[%(facet)s][1] - x%(restriction)s[v0][1])  + a2*(x%(restriction)s[%(facet)s][2] - x%(restriction)s[v0][2]) < 0;"""

_facet_normal_1D = """
// Facet normals are 1.0 or -1.0:   (-1.0) <-- X------X --> (1.0)
const double n%(restriction)s = %(direction)sdirection ? 1.0 : -1.0;"""

_facet_normal_2D = """\
// Compute facet normals from the facet scale factor constants
const double n%(restriction)s0 = %(direction)sdirection ? dx1 / det : -dx1 / det;
const double n%(restriction)s1 = %(direction)sdirection ? -dx0 / det : dx0 / det;"""

_facet_normal_3D = """\
// Compute facet normals from the facet scale factor constants
const double n%(restriction)s0 = %(direction)sdirection ? a0 / det : -a0 / det;
const double n%(restriction)s1 = %(direction)sdirection ? a1 / det : -a1 / det;
const double n%(restriction)s2 = %(direction)sdirection ? a2 / det : -a2 / det;"""

_cell_volume_1D = """\
// Cell Volume.
const double volume%(restriction)s = std::abs(detJ%(restriction)s);"""

_cell_volume_2D = """\
// Cell Volume.
const double volume%(restriction)s = std::abs(detJ%(restriction)s)/2.0;"""

_cell_volume_3D = """\
// Cell Volume.
const double volume%(restriction)s = std::abs(detJ%(restriction)s)/6.0;"""

_circumradius_1D = """\
// Compute circumradius, in 1D it is equal to the cell volume.
const double circumradius%(restriction)s = std::abs(detJ%(restriction)s);"""

_circumradius_2D = """\
// Compute circumradius, assuming triangle is embedded in 2D.
const double v1v2%(restriction)s  = std::sqrt( (x%(restriction)s[2][0] - x%(restriction)s[1][0])*(x%(restriction)s[2][0] - x%(restriction)s[1][0]) + (x%(restriction)s[2][1] - x%(restriction)s[1][1])*(x%(restriction)s[2][1] - x%(restriction)s[1][1]) );
const double v0v2%(restriction)s  = std::sqrt( J%(restriction)s_11*J%(restriction)s_11 + J%(restriction)s_01*J%(restriction)s_01 );
const double v0v1%(restriction)s  = std::sqrt( J%(restriction)s_00*J%(restriction)s_00 + J%(restriction)s_10*J%(restriction)s_10 );

const double circumradius%(restriction)s = 0.25*(v1v2%(restriction)s*v0v2%(restriction)s*v0v1%(restriction)s)/(volume%(restriction)s);"""

_circumradius_3D = """\
// Compute circumradius.
const double v1v2%(restriction)s  = std::sqrt( (x%(restriction)s[2][0] - x%(restriction)s[1][0])*(x%(restriction)s[2][0] - x%(restriction)s[1][0]) + (x%(restriction)s[2][1] - x%(restriction)s[1][1])*(x%(restriction)s[2][1] - x%(restriction)s[1][1]) + (x%(restriction)s[2][2] - x%(restriction)s[1][2])*(x%(restriction)s[2][2] - x%(restriction)s[1][2]) );
const double v0v2%(restriction)s  = std::sqrt(J%(restriction)s_01*J%(restriction)s_01 + J%(restriction)s_11*J%(restriction)s_11 + J%(restriction)s_21*J%(restriction)s_21);
const double v0v1%(restriction)s  = std::sqrt(J%(restriction)s_00*J%(restriction)s_00 + J%(restriction)s_10*J%(restriction)s_10 + J%(restriction)s_20*J%(restriction)s_20);
const double v0v3%(restriction)s  = std::sqrt(J%(restriction)s_02*J%(restriction)s_02 + J%(restriction)s_12*J%(restriction)s_12 + J%(restriction)s_22*J%(restriction)s_22);
const double v1v3%(restriction)s  = std::sqrt( (x%(restriction)s[3][0] - x%(restriction)s[1][0])*(x%(restriction)s[3][0] - x%(restriction)s[1][0]) + (x%(restriction)s[3][1] - x%(restriction)s[1][1])*(x%(restriction)s[3][1] - x%(restriction)s[1][1]) + (x%(restriction)s[3][2] - x%(restriction)s[1][2])*(x%(restriction)s[3][2] - x%(restriction)s[1][2]) );
const double v2v3%(restriction)s  = std::sqrt( (x%(restriction)s[3][0] - x%(restriction)s[2][0])*(x%(restriction)s[3][0] - x%(restriction)s[2][0]) + (x%(restriction)s[3][1] - x%(restriction)s[2][1])*(x%(restriction)s[3][1] - x%(restriction)s[2][1]) + (x%(restriction)s[3][2] - x%(restriction)s[2][2])*(x%(restriction)s[3][2] - x%(restriction)s[2][2]) );
const  double la%(restriction)s   = v1v2%(restriction)s*v0v3%(restriction)s;
const  double lb%(restriction)s   = v0v2%(restriction)s*v1v3%(restriction)s;
const  double lc%(restriction)s   = v0v1%(restriction)s*v2v3%(restriction)s;
const  double s%(restriction)s    = 0.5*(la%(restriction)s+lb%(restriction)s+lc%(restriction)s);
const  double area%(restriction)s = std::sqrt(s%(restriction)s*(s%(restriction)s-la%(restriction)s)*(s%(restriction)s-lb%(restriction)s)*(s%(restriction)s-lc%(restriction)s));

const double circumradius%(restriction)s = area%(restriction)s / ( 6.0*volume%(restriction)s );"""

_facet_area_1D = """\
// Facet Area (FIXME: Should this be 0.0?).
const double facet_area = 1.0;"""

_facet_area_2D = """\
// Facet Area.
const double facet_area = det;"""

_facet_area_3D = """\
// Facet Area (divide by two because 'det' is scaled by area of reference triangle).
const double facet_area = det/2.0;"""

evaluate_basis_dofmap = """\
unsigned int element = 0;
unsigned int tmp = 0;
for (unsigned int j = 0; j < %d; j++)
{
  if (tmp +  dofs_per_element[j] > i)
  {
    i -= tmp;
    element = element_types[j];
    break;
  }
  else
    tmp += dofs_per_element[j];
}"""

# Used in evaluate_basis_derivatives. For second order derivatives in 2D it will
# generate the combinations: [(0, 0), (0, 1), (1, 0), (1, 1)] (i.e., xx, xy, yx, yy)
# which will also be the ordering of derivatives in the return value.
combinations_snippet = """\
// Declare pointer to two dimensional array that holds combinations of derivatives and initialise
unsigned int **%(combinations)s = new unsigned int *[%(num_derivatives)s];
for (unsigned int row = 0; row < %(num_derivatives)s; row++)
{
  %(combinations)s[row] = new unsigned int [%(n)s];
  for (unsigned int col = 0; col < %(n)s; col++)
    %(combinations)s[row][col] = 0;
}

// Generate combinations of derivatives
for (unsigned int row = 1; row < %(num_derivatives)s; row++)
{
  for (unsigned int num = 0; num < row; num++)
  {
    for (unsigned int col = %(n)s-1; col+1 > 0; col--)
    {
      if (%(combinations)s[row][col] + 1 > %(topological_dimension-1)s)
        %(combinations)s[row][col] = 0;
      else
      {
        %(combinations)s[row][col] += 1;
        break;
      }
    }
  }
}"""

_transform_interval_snippet = """\
// Compute inverse of Jacobian
const double %(K)s[1][1] =  {{K_00}};

// Declare transformation matrix
// Declare pointer to two dimensional array and initialise
double **%(transform)s = new double *[%(num_derivatives)s];

for (unsigned int j = 0; j < %(num_derivatives)s; j++)
{
  %(transform)s[j] = new double [%(num_derivatives)s];
  for (unsigned int k = 0; k < %(num_derivatives)s; k++)
    %(transform)s[j][k] = 1;
}

// Construct transformation matrix
for (unsigned int row = 0; row < %(num_derivatives)s; row++)
{
  for (unsigned int col = 0; col < %(num_derivatives)s; col++)
  {
    for (unsigned int k = 0; k < %(n)s; k++)
      %(transform)s[row][col] *= %(K)s[%(combinations)s[col][k]][%(combinations)s[row][k]];
  }
}"""

_transform_triangle_snippet = """\
// Compute inverse of Jacobian
const double %(K)s[2][2] = \
{{K_00, K_01},\
 {K_10, K_11}};

// Declare transformation matrix
// Declare pointer to two dimensional array and initialise
double **%(transform)s = new double *[%(num_derivatives)s];

for (unsigned int j = 0; j < %(num_derivatives)s; j++)
{
  %(transform)s[j] = new double [%(num_derivatives)s];
  for (unsigned int k = 0; k < %(num_derivatives)s; k++)
    %(transform)s[j][k] = 1;
}

// Construct transformation matrix
for (unsigned int row = 0; row < %(num_derivatives)s; row++)
{
  for (unsigned int col = 0; col < %(num_derivatives)s; col++)
  {
    for (unsigned int k = 0; k < %(n)s; k++)
      %(transform)s[row][col] *= %(K)s[%(combinations)s[col][k]][%(combinations)s[row][k]];
  }
}"""

_transform_tetrahedron_snippet = """\
// Compute inverse of Jacobian
const double %(K)s[3][3] = \
{{K_00, K_01, K_02},\
 {K_10, K_11, K_12},\
 {K_20, K_21, K_22}};

// Declare transformation matrix
// Declare pointer to two dimensional array and initialise
double **%(transform)s = new double *[%(num_derivatives)s];

for (unsigned int j = 0; j < %(num_derivatives)s; j++)
{
  %(transform)s[j] = new double [%(num_derivatives)s];
  for (unsigned int k = 0; k < %(num_derivatives)s; k++)
    %(transform)s[j][k] = 1;
}

// Construct transformation matrix
for (unsigned int row = 0; row < %(num_derivatives)s; row++)
{
  for (unsigned int col = 0; col < %(num_derivatives)s; col++)
  {
    for (unsigned int k = 0; k < %(n)s; k++)
      %(transform)s[row][col] *= %(K)s[%(combinations)s[col][k]][%(combinations)s[row][k]];
  }
}"""

# Codesnippets used in evaluate_dof
_map_onto_physical_1D = """\
// Evaluate basis functions for affine mapping
const double w0 = 1.0 - X_%(i)d[%(j)s][0];
const double w1 = X_%(i)d[%(j)s][0];

// Compute affine mapping y = F(X)
y[0] = w0*x[0][0] + w1*x[1][0];"""

_map_onto_physical_2D = """\
// Evaluate basis functions for affine mapping
const double w0 = 1.0 - X_%(i)d[%(j)s][0] - X_%(i)d[%(j)s][1];
const double w1 = X_%(i)d[%(j)s][0];
const double w2 = X_%(i)d[%(j)s][1];

// Compute affine mapping y = F(X)
y[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0];
y[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1];"""

_map_onto_physical_3D = """\
// Evaluate basis functions for affine mapping
const double w0 = 1.0 - X_%(i)d[%(j)s][0] - X_%(i)d[%(j)s][1] - X_%(i)d[%(j)s][2];
const double w1 = X_%(i)d[%(j)s][0];
const double w2 = X_%(i)d[%(j)s][1];
const double w3 = X_%(i)d[%(j)s][2];

// Compute affine mapping y = F(X)
y[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0] + w3*x[3][0];
y[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1] + w3*x[3][1];
y[2] = w0*x[0][2] + w1*x[1][2] + w2*x[2][2] + w3*x[3][2];"""

_ip_coordinates_1D = """\
X%(num_ip)d[0] = %(name)s[%(ip)s][0]*x%(restriction)s[0][0] + \
%(name)s[%(ip)s][1]*x%(restriction)s[1][0];"""

_ip_coordinates_2D = """\
X%(num_ip)d[0] = %(name)s[%(ip)s][0]*x%(restriction)s[0][0] + \
%(name)s[%(ip)s][1]*x%(restriction)s[1][0] + %(name)s[%(ip)s][2]*x%(restriction)s[2][0];
X%(num_ip)d[1] = %(name)s[%(ip)s][0]*x%(restriction)s[0][1] + \
%(name)s[%(ip)s][1]*x%(restriction)s[1][1] + %(name)s[%(ip)s][2]*x%(restriction)s[2][1];"""

_ip_coordinates_3D = """\
X%(num_ip)d[0] = %(name)s[%(ip)s][0]*x%(restriction)s[0][0] + \
%(name)s[%(ip)s][1]*x%(restriction)s[1][0] + %(name)s[%(ip)s][2]*x%(restriction)s[2][0] + \
%(name)s[%(ip)s][3]*x%(restriction)s[3][0];
X%(num_ip)d[1] = %(name)s[%(ip)s][0]*x%(restriction)s[0][1] + \
%(name)s[%(ip)s][1]*x%(restriction)s[1][1] + %(name)s[%(ip)s][2]*x%(restriction)s[2][1] + \
%(name)s[%(ip)s][3]*x%(restriction)s[3][1];
X%(num_ip)d[2] = %(name)s[%(ip)s][0]*x%(restriction)s[0][2] + \
%(name)s[%(ip)s][1]*x%(restriction)s[1][2] + %(name)s[%(ip)s][2]*x%(restriction)s[2][2] + \
%(name)s[%(ip)s][3]*x%(restriction)s[3][2];"""

# Codesnippets used in evaluatebasis[|derivatives]
_map_coordinates_FIAT_interval = """\
// Get coordinates and map to the reference (FIAT) element
double X = (2.0*coordinates[0] - x[0][0] - x[1][0]) / J_00;"""

_map_coordinates_FIAT_interval_in_2D = """\
// Get coordinates and map to the reference (FIAT) element
double X = 2*(std::sqrt(std::pow(coordinates[0]-x[0][0], 2) + 
                        std::pow(coordinates[1]-x[0][1], 2))/ detJ) - 1.0;"""

_map_coordinates_FIAT_triangle = """\
// Compute constants
const double C0 = x[1][0] + x[2][0];
const double C1 = x[1][1] + x[2][1];

// Get coordinates and map to the reference (FIAT) element
double X = (J_01*(C1 - 2.0*coordinates[1]) + J_11*(2.0*coordinates[0] - C0)) / detJ;
double Y = (J_00*(2.0*coordinates[1] - C1) + J_10*(C0 - 2.0*coordinates[0])) / detJ;"""

_map_coordinates_FIAT_tetrahedron = """\
// Compute constants
const double C0 = x[3][0] + x[2][0] + x[1][0] - x[0][0];
const double C1 = x[3][1] + x[2][1] + x[1][1] - x[0][1];
const double C2 = x[3][2] + x[2][2] + x[1][2] - x[0][2];

// Get coordinates and map to the reference (FIAT) element
double X = (d_00*(2.0*coordinates[0] - C0) + d_10*(2.0*coordinates[1] - C1) + d_20*(2.0*coordinates[2] - C2)) / detJ;
double Y = (d_01*(2.0*coordinates[0] - C0) + d_11*(2.0*coordinates[1] - C1) + d_21*(2.0*coordinates[2] - C2)) / detJ;
double Z = (d_02*(2.0*coordinates[0] - C0) + d_12*(2.0*coordinates[1] - C1) + d_22*(2.0*coordinates[2] - C2)) / detJ;
"""


# Mappings to code snippets used by format

jacobian = {1: {1:_jacobian_1D}, 
            2: {2:_jacobian_2D, 1:_jacobian_2D_1D},
            3: {3:_jacobian_3D, 2:_jacobian_3D_2D, 1:_jacobian_3D_1D}}

inverse_jacobian = {1: {1:_inverse_jacobian_1D},
                    2: {2:_inverse_jacobian_2D, 1:_inverse_jacobian_2D_1D},
                    3: {3:_inverse_jacobian_3D}}

facet_determinant = {1: _facet_determinant_1D,
                     2: _facet_determinant_2D,
                     3: _facet_determinant_3D}

map_onto_physical = {1: _map_onto_physical_1D,
                     2: _map_onto_physical_2D,
                     3: _map_onto_physical_3D}

fiat_coordinate_map = {"interval": {1:_map_coordinates_FIAT_interval,
                                    2:_map_coordinates_FIAT_interval_in_2D},
                       "triangle": {2:_map_coordinates_FIAT_triangle},
                       "tetrahedron": {3:_map_coordinates_FIAT_tetrahedron}}

transform_snippet = {"interval": _transform_interval_snippet,
                     "triangle": _transform_triangle_snippet,
                     "tetrahedron": _transform_tetrahedron_snippet}

normal_direction = {1: _normal_direction_1D,
                    2: _normal_direction_2D,
                    3: _normal_direction_3D}

facet_normal = {1: _facet_normal_1D,
                2: _facet_normal_2D,
                3: _facet_normal_3D}

ip_coordinates = {1: (3, _ip_coordinates_1D),
                  2: (10, _ip_coordinates_2D),
                  3: (21, _ip_coordinates_3D)}

cell_volume = {1: _cell_volume_1D,
               2: _cell_volume_2D,
               3: _cell_volume_3D}

circumradius = {1: _circumradius_1D,
                2: _circumradius_2D,
                3: _circumradius_3D}

facet_area = {1: _facet_area_1D,
              2: _facet_area_2D,
              3: _facet_area_3D}

