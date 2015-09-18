"This module provides simple C++ code for verification of UFC code."

# Copyright (C) 2010 Kristian B. Oelgaard
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
# First added:  2010-01-18
# Last changed: 2010-01-18

#evaluate_basis_code = """\
##include <iostream>
##include <ufc.h>
##include "test.h"

#int main()
#{
#  // Create element
#  %(element)s element;

#  // Size of dof_values
#  // FIXME: This will not work for TensorElements
#  int N = element.value_dimension(0);

#  // Create values
#  double* dof_values = new double[N];
#  for (unsigned int i = 0; i < N; i++)
#    dof_values[i] = 0.0;

#  // Create vertex coordinates and fill with some arbitrary data
#  double coordinate_dofs[24] = {0.90, 0.34, 0.45,
#                                   0.56, 0.76, 0.83,
#                                   0.98, 0.78, 0.19,
#                                   0.12, 0.56, 0.66,
#                                   0.96, 0.78, 0.63,
#                                   0.11, 0.35, 0.49,
#                                   0.51, 0.88, 0.65,
#                                   0.98, 0.45, 0.01};

#  // Random coordinates where we want to evaluate the basis functions
#  double coordinates[3] = {0.32, 0.51, 0.05};

#  // Loop element space dimension and call evaluate_basis.
#  for (unsigned int i = 0; i < element.space_dimension(); i++)
#  {
#    element.evaluate_basis(i, dof_values, coordinates, coordinate_dofs, 0);
#    // Print values
#    for (unsigned int j = 0; j < N; j++)
#      std::cout << dof_values[j] << " ";
#  }
#  std::cout << std::endl;
#  return 0;
#}
#"""

evaluate_basis_code_fiat = """\
#include <iostream>
#include <ufc.h>
#include <cstdlib>
#include "test.h"

int main(int argc, char* argv[])
{
  // Create element
  %(element)s element;


  // Get derivative order
  unsigned int n = std::atoi(argv[1]);

  // Value dimension
  int N;
  if (element.value_rank() == 0)
  {
    N = 1;
  }
  else
  {
    N = 1;
    for (unsigned int i = 0; i < element.value_rank(); i++)
    {
      N = N * element.value_dimension(i);
    }
  }

  // Compute number of derivatives.
  unsigned int  num_derivatives = 1;
  for (unsigned int r = 0; r < n; r++)
  {
    num_derivatives *= %(dim)d;
  }

  // Create values
  unsigned int num_dof_vals = N*num_derivatives;
  double* dof_values = new double[num_dof_vals];
  for (unsigned int i = 0; i < num_dof_vals; i++)
    dof_values[i] = 0.0;

  %(cell_ref_coords)s
  double coordinate_dofs[%(num_coords)d*%(dim)d];
  int k = 0;
  for (int i = 0; i < %(num_coords)d; i++)
  {
    for (int j = 0; j < %(dim)d; j++)
      coordinate_dofs[k++] = cell_ref_coords[i][j];
  }

  // Random points where we want to evaluate the basis functions
  // coordinates of dofs and three arbitrary points on the reference cell.
  double points[%(num_points)d][%(dim)d] = %(points)s

  // Init array of coordinates
  double coordinates[3] = {0,0,0};

  std::cout.precision(8);
  std::cout.setf(std::ios::fixed);

  // If we're testing evaluate_basis, loop all points.
  if (n == 0)
  {
    for (unsigned int p = 0; p < %(num_points)d; p++)
    {
      for (unsigned int d = 0; d < %(dim)d; d++)
      {
        coordinates[d] = points[p][d];
      }
      // Loop element space dimension and call evaluate_basis.
      for (unsigned int i = 0; i < element.space_dimension(); i++)
      {
        element.evaluate_basis(i, dof_values, coordinates, coordinate_dofs, 0);

        // Print values
        for (unsigned int j = 0; j < num_dof_vals; j++)
          std::cout << dof_values[j] << " ";
      }
      std::cout << std::endl;
    }
  }
  else
  {
    // Else loop the arbitrary 3 points, otherwise the number of tests explode
    // with the element.space_dimension()^2.
    for (unsigned int p = element.space_dimension(); p < %(num_points)d; p++)
    {
      for (unsigned int d = 0; d < %(dim)d; d++)
      {
        coordinates[d] = points[p][d];
      }
      // Loop element space dimension and call evaluate_basis.
      for (unsigned int i = 0; i < element.space_dimension(); i++)
      {
        element.evaluate_basis_derivatives(i, n, dof_values, coordinates, coordinate_dofs, 0);

        // Print values
        for (unsigned int j = 0; j < num_dof_vals; j++)
          std::cout << dof_values[j] << " ";
      }
      std::cout << std::endl;
    }
  }

  return 0;
}
"""
