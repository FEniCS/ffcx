// Copyright (C) 2003-2009 Anders Logg
//
// This file is part of FFC.
//
// DOLFIN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DOLFIN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN.  If not, see <http://www.gnu.org/licenses/>.
//
// First added:  2012-08-20
// Last changed: 2012-09-05

#ifndef __LEGENDRE_H
#define __LEGENDRE_H

/// Legendre polynomial of given degree n on the interval [-1,1].
///
///   P0(x) = 1
///   P1(x) = x
///   P2(x) = (3x^2 - 1) / 2
///   ...
///
/// The function values and derivatives are computed using
/// three-term recurrence formulas.

#include <vector>
#define EPSILON 10e-15

class Legendre
{
 public:

  Legendre(unsigned int n);

  /// Evaluation at given point
  double operator() (double x);

  /// Evaluation of derivative at given point
  double ddx(double x);

  /// Evaluation of second derivative at given point
  double d2dx(double x);

  /// Evaluation of arbitrary order, nn <= n (useful ie in RadauQuadrature)
  double eval(unsigned int nn, double x);

  double ddx(unsigned int n, double x);
  double d2dx(unsigned int n, double x);
 
 private:

  const unsigned int n;
  double cache_x;
  std::vector<double> cache;

};
#endif
