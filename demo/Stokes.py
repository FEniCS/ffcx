# Copyright (C) 2005-2007 Anders Logg
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
# Modified by: Martin Sandve Alnes (2009)
# Date: 2009-03-02
#
# The bilinear form a(v, u) and Linear form L(v) for the Stokes
# equations using a mixed formulation (Taylor-Hood elements).
from ufl import (Coefficient, FiniteElement, TestFunctions, TrialFunctions,
                 VectorElement, div, dot, dx, grad, inner, triangle)

cell = triangle
P2 = VectorElement("Lagrange", cell, 2)
P1 = FiniteElement("Lagrange", cell, 1)
TH = P2 * P1

(v, q) = TestFunctions(TH)
(u, p) = TrialFunctions(TH)

f = Coefficient(P2)

a = (inner(grad(v), grad(u)) - div(v) * p + q * div(u)) * dx
L = dot(v, f) * dx
