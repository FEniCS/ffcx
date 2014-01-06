// Copyright (C) 2003-2008 Anders Logg
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
// Modified by Benjamin Kehlet 2011-2012
//
// First added:  2012-08-20
// Last changed: 2012-09-05


#include "Legendre.h"
#include <cmath>

//-----------------------------------------------------------------------------
Legendre::Legendre(unsigned int n) : n(n), cache_x(0.0), cache(n + 1)
{
  cache[0] = 1.0; //constant value

  // eval to initialize cache
  eval(n, -1.0);
}
//-----------------------------------------------------------------------------
double Legendre::operator() (double x)
{
  return eval(n, x);
}
//-----------------------------------------------------------------------------
double Legendre::ddx(double x)
{
  return ddx(n, x);
}
//-----------------------------------------------------------------------------
double Legendre::d2dx(double x)
{
  return d2dx(n, x);
}
//-----------------------------------------------------------------------------
double Legendre::eval(unsigned int nn, double x)
{
  //recursive formula, BETA page 254
  //return ( (2.0*nn-1.0)*x*eval(nn-1, x) - (nn-1.0)*eval(nn-2, x) ) / nn;

  //The special cases
  if (n == 0)
    return 1.0;
  else if (n == 1)
    return x;

  //check cache
  if (x != cache_x)
  {
    cache[1] = x;
    for (unsigned int i = 2; i <= n; ++i)
    {
      double ii = i;
      cache[i] = ( (2.0*ii-1.0)*x*cache[i-1] - (ii-1.0)*cache[i-2] ) / ii;
    }
    cache_x = x;
  }

  return cache[nn];
}
//-----------------------------------------------------------------------------
double Legendre::ddx(unsigned int n, double x)
{
  // Special cases
  if (n == 0)
    return 0.0;
  else if (n == 1)
    return 1.0;

  // Avoid division by zero
  if (std::abs(x - 1.0) < EPSILON)
    x -= 2.0*EPSILON;

  if (std::abs(x + 1.0) < EPSILON)
    x += 2.0*EPSILON;

  // Formula, BETA page 254
  const double nn = n;
  return nn * (x*eval(n, x) - eval(n-1, x)) / (x*x - 1.0);
}
//-----------------------------------------------------------------------------
double Legendre::d2dx(unsigned int, double x)
{
  // Special case n = 0
  if (n == 0)
    return 0.0;

  // Special case n = 1
  if (n == 1)
    return 0.0;

  // Avoid division by zero
  if (std::abs(x - 1.0) < EPSILON)
    x -= 2.0*EPSILON;
  if (std::abs(x + 1.0) < EPSILON)
    x += 2.0*EPSILON;

  // Formula, BETA page 254
  const double nn = double(n);
  return (2.0*x*ddx(n, x) - nn*(nn+1)*eval(n, x)) / (1.0-x*x);
}
//-----------------------------------------------------------------------------
