"This module provides simple C++ code for verification of UFC code."

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@gmail.com)"
__date__ = "2010-01-18"
__copyright__ = "Copyright (C) 2010 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Common code for all integral types
evaluate_basis_code = """\
#include <iostream>
#include <ufc.h>
#include "test.h"

int main()
{
  // Create element
  %(element)s element;

  // Size of dof_values
  // FIXME: This will not work for TensorElements
  int N = element.value_dimension(0);

  // Create values
  double* dof_values = new double[N];
  for (int i = 0; i < N; i++)
    dof_values[i] = 0.0;

  // Create cell and fill with some arbitrary data
  double cell_coordinates[8][3] = {{0.90, 0.34, 0.45},
                                   {0.56, 0.76, 0.83},
                                   {0.98, 0.78, 0.19},
                                   {0.12, 0.56, 0.66},
                                   {0.96, 0.78, 0.63},
                                   {0.11, 0.35, 0.49},
                                   {0.51, 0.88, 0.65},
                                   {0.98, 0.45, 0.01}};
  ufc::cell cell;
  cell.coordinates = new double * [8];
  for (int i = 0; i < 8; i++)
  {
    cell.coordinates[i] = new double[3];
    for (int j = 0; j < 3; j++)
      cell.coordinates[i][j] = cell_coordinates[i][j];
  }

  // Random coordinates where we want to evaluate the basis functions
  double coordinates[3] = {0.32, 0.51, 0.05};

  // Loop element space dimension and call evaluate_basis.
  for (unsigned int i = 0; i < element.space_dimension(); i++)
  {
    element.evaluate_basis(i, dof_values, coordinates, cell);
    // Print values
    for (int j = 0; j < N; j++)
      std::cout << dof_values[j] << " ";
  }
  std::cout << std::endl;
  return 0;
}
"""

