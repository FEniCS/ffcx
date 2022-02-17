# Copyright (C) 2008 Kristian B. Oelgaard
#
# This file is part of UFL.
#
# UFL is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# UFL is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with UFL. If not, see <http://www.gnu.org/licenses/>.
#
# First added:  2008-03-31
# Last changed: 2008-03-31
#
# The linearised bilinear form a(u,v) and linear form L(v) for
# the nonlinear equation - div (1+u) grad u = f (non-linear Poisson)
from ufl import (Coefficient, FiniteElement, TestFunction, TrialFunction,
                 VectorElement, dot, dx, grad, i, triangle)

element = FiniteElement("Lagrange", triangle, 2)

QE = FiniteElement("Quadrature", triangle, 2, quad_scheme="default")
sig = VectorElement("Quadrature", triangle, 1, quad_scheme="default")

v = TestFunction(element)
u = TrialFunction(element)
u0 = Coefficient(element)
C = Coefficient(QE)
sig0 = Coefficient(sig)
f = Coefficient(element)

a = v.dx(i) * C * u.dx(i) * dx(metadata={"quadrature_degree": 2}) + v.dx(i) * 2 * u0 * u * u0.dx(i) * dx
L = v * f * dx - dot(grad(v), sig0) * dx(metadata={"quadrature_degree": 1})
