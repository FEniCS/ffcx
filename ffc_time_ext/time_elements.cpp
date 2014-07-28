// Copyright (C) 2012 Benjamin Kehlet
//
// This file is part of FFC.
//
// FFC is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// FFC is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with FFC.  If not, see <http://www.gnu.org/licenses/>.
//
// First added:  2012-08-20
// Last changed: 2012-09-05

#include <iostream>
#include "LobattoQuadrature.h"
#include "RadauQuadrature.h"
#include "Legendre.h"

void compute_lobatto_points(double* points, const unsigned int degree)
{
  // Compute the nodal basis
  LobattoQuadrature lobatto(degree + 1);
  for (unsigned int i = 0; i < degree +1; i++)
    points[i] = (lobatto.points[i] + 1.0) / 2.0;
}

void compute_radau_points  (double* points, const unsigned int degree)
{
  RadauQuadrature radau(degree+1);

  for (unsigned int i = 0; i < degree+1; i++)
    points[degree-i] = (-radau.points[i] + 1.0) / 2.0;
}

void compute_legendre_coeffs(double* coeffs, const double *points, const unsigned int num_points, const unsigned int degree)
{
  for (unsigned int i = 0; i < degree; i++)
  {
    Legendre p(i);
    for (unsigned int j = 0; j < num_points; j++)
    {
      coeffs[i*num_points + j] = p(points[j]*2.0 -1.0);
    }
  }
}
