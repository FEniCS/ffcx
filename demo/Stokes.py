# Copyright (C) 2005-2007 Anders Logg
#
# This file is part of FFCx.
#
# FFCx is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFCx is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFCx. If not, see <http://www.gnu.org/licenses/>.
#
# The bilinear form a(u, v) and Linear form L(v) for the Stokes
# equations using a mixed formulation (Taylor-Hood elements).

# Compile this form with FFCx: ffcx Stokes.ufl
from ufl import (Coefficient, FiniteElement, TestFunctions, TrialFunctions,
                 VectorElement, div, dx, grad, inner, triangle)

P2 = VectorElement("Lagrange", triangle, 2)
P1 = FiniteElement("Lagrange", triangle, 1)
TH = P2 * P1

(u, p) = TrialFunctions(TH)
(v, q) = TestFunctions(TH)

f = Coefficient(P2)

a = (inner(grad(u), grad(v)) - div(v) * p + div(u) * q) * dx
L = inner(f, v) * dx
