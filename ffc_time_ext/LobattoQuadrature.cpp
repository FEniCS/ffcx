// Copyright (C) 2003-2006 Anders Logg
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

#include "LobattoQuadrature.h"
#include "Legendre.h"

#include <cmath>

//-----------------------------------------------------------------------------
LobattoQuadrature::LobattoQuadrature(unsigned int n)
  : points(n)
{
  // FIXME: Do proper arguement checking
  // if (n < 2)
  //   error("Lobatto quadrature requires at least 2 points.");

  // init();

  // if (!check(2*n - 3))
  //   error("Lobatto quadrature not ok, check failed.");



  // Compute the Lobatto quadrature points in [-1,1] as the endpoints
  // and the zeroes of the derivatives of the Legendre polynomials
  // using Newton's method

  //const unsigned int n = points.size();

  // Special case n = 1 (should not be used)
  if (n == 1)
  {
    points[0] = 0.0;
    return;
  }

  // Special case n = 2
  if (n == 2)
  {
    points[0] = -1.0;
    points[1] = 1.0;
    return;
  }

  Legendre p(n - 1);
  double x, dx;

  // Set the first and last nodal points which are 0 and 1
  points[0] = -1.0;
  points[n - 1] = 1.0;

  // Compute the rest of the nodes by Newton's method
  for (unsigned int i = 1; i <= ((n-1)/2); i++)
  {

    // Initial guess
    x = cos(3.1415926*double(i)/double(n - 1));

    // Newton's method
    do
    {
      dx = -p.ddx(x)/p.d2dx(x);
      x  = x + dx;
    } while (std::abs(dx) > EPSILON);

    // Save the value using the symmetry of the points
    points[i] = -x;
    points[n - 1 - i] = x;

  }

  // Fix the middle node
  if ((n % 2) != 0)
    points[n/2] = 0.0;
}

