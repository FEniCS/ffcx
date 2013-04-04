// Copyright (C) 2003-2009 Anders Logg
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

#ifndef __RADAU_QUADRATURE_H
#define __RADAU_QUADRATURE_H

/// Radau (Gauss-Radau) quadrature on the interval [-1,1].
/// The n quadrature points are given by the zeros of
///
///     ( Pn-1(x) + Pn(x) ) / (1+x)
///
/// where Pn is the n:th Legendre polynomial.
///
/// The quadrature points are computed using Newton's method, and
/// the quadrature weights are computed by solving a linear system
/// determined by the condition that Radau quadrature with n points
/// should be exact for polynomials of degree 2n-2.

#include <vector>

class RadauQuadrature
{
 public:

  /// Create Radau quadrature with n points
  RadauQuadrature(unsigned int n);

  std::vector<double> points;

};

#endif
