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

#include "RadauQuadrature.h"
#include "Legendre.h"

#include <cmath>

//-----------------------------------------------------------------------------
RadauQuadrature::RadauQuadrature(unsigned int n)
  : points(n+1)
{
  // Compute the Radau quadrature points in [-1,1] as -1 and the zeros
  // of ( Pn-1(x) + Pn(x) ) / (1+x) where Pn is the n:th Legendre
  // polynomial. Computation is a little different than for Gauss and
  // Lobatto quadrature, since we don't know of any good initial
  // approximation for the Newton iterations.

  // Special case n = 1
  if (n == 1)
  {
    points[0] = -1.0;
    return;
  }

  Legendre p(n);
  double x, dx, step, sign;

  // Set size of stepping for seeking starting points
  step = 1.0/(double(n - 1)*15.0);

  // Set the first nodal point which is -1
  points[0] = -1.0;

  // Start at -1 + step
  x = -1.0 + step;

  // Set the sign at -1 + epsilon
  sign = ((p.eval(n - 1, x) + p(x)) > 0 ? 1.0 : -1.0);

  // Compute the rest of the nodes by Newton's method
  for (unsigned int i = 1; i < n; i++)
  {

    // Step to a sign change
    while ((p.eval(n - 1, x) + p(x))*sign > 0.0)
      x += step;

    // Newton's method
    do
    {
      dx = -(p.eval(n-1, x) + p(x))/(p.ddx(n - 1, x) + p.ddx(x));
      x  = x + dx;
    } while (std::abs(dx) > EPSILON);

    // Set the node value
    points[i] = x;

    // Fix step so that it's not too large
    if (step > (points[i] - points[i-1])/10.0)
      step = (points[i] - points[i-1])/10.0;

    // Step forward
    sign = -sign;
    x += step;
  }
}
