"This module provides simple C++ code for verification of UFC code."

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2009-03-15 -- 2009-03-15"
__copyright__ = "Copyright (C) 2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Integral types
integral_types = ["cell_integral", "exterior_facet_integral", "interior_facet_integral"]

# Common code for all integral types
tabulate_tensor_code_common = """\
#include <iostream>
#include <ufc.h>
#include "%(header)s"

int main()
{
  // Size of arrays
  int n = %(n)s;
  int N = %(N)s;

  // Create integral
  %(integral)s integral;

  // Create tensor
  double* A = new double[N];
  for (int i = 0; i < N; i++)
    A[i] = 0.0;

  // Create coefficients and fill with some arbitrary data
  double** w = new double * [n];
  for (int i = 0; i < n; i++)
  {
    w[i] = new double[N];
    for (int j = 0; j < N; j++)
      w[i][j] = 1.0 / static_cast<double>(i + j + 1);
  }

  // Create cell and fill with some arbitrary data
  double coordinates[8][3] = {{0.90, 0.34, 0.45},
                              {0.56, 0.76, 0.83},
                              {0.98, 0.78, 0.19},
                              {0.12, 0.56, 0.66},
                              {0.96, 0.78, 0.63},
                              {0.11, 0.35, 0.49},
                              {0.51, 0.88, 0.65},
                              {0.98, 0.45, 0.01}};
  ufc::cell cell0;
  ufc::cell cell1;
  cell0.coordinates = new double * [8];
  cell1.coordinates = new double * [8];
  for (int i = 0; i < 8; i++)
  {
    cell0.coordinates[i] = new double[3];
    cell1.coordinates[i] = new double[3];
    for (int j = 0; j < 3; j++)
    {
      cell0.coordinates[i][j] = coordinates[i][j];
      cell1.coordinates[i][j] = coordinates[i][j] + 0.1;
    }
  }

  // Call tabulate_tensor
  <tabulate_tensor>;
  
  // Print values
  for (int i = 0; i < N; i++)
    std::cout << A[i] << " ";

  return 0;
}
"""

# Code for different integral types
i = integral_types
c = tabulate_tensor_code_common
tabulate_tensor_code = {i[0] : c.replace("<tabulate_tensor>", "integral.tabulate_tensor(A, w, cell0)"),
                        i[1] : c.replace("<tabulate_tensor>", "integral.tabulate_tensor(A, w, cell0, 0)"),
                        i[2] : c.replace("<tabulate_tensor>", "integral.tabulate_tensor(A, w, cell0, cell1, 0, 0)")}
