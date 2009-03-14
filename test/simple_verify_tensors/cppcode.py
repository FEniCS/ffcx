"This module provides simple C++ code for verification of UFC code."

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2009-03-15 -- 2009-03-15"
__copyright__ = "Copyright (C) 2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

tabulate_tensor_code = """\
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
  ufc::cell cell;
  cell.coordinates = new double * [8];
  for (int i = 0; i < 8; i++)
  {
    cell.coordinates[i] = new double[3];
    for (int j = 0; j < 3; j++)
      cell.coordinates[i][j] = coordinates[i][j];
  }

  // Call tabulate_tensor
  integral.tabulate_tensor(A, w, cell);
  
  // Print values
  for (int i = 0; i < N; i++)
    std::cout << A[i] << " ";

  return 0;
}
"""
