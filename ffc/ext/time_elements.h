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


void compute_lobatto_points(double* points, const unsigned int degree);

void compute_radau_points  (double* points, const unsigned int degree);

void compute_legendre_coeffs(double* coeffs, 
			     const double *points, 
			     const unsigned int num_points, 
			     const unsigned int degree);
