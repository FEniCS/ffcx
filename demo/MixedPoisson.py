# Copyright (C) 2006-2009 Anders Logg and Marie Rognes
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
# Modified by Martin Sandve Alnes, 2009
#
# Last changed: 2009-01-12
#
# The bilinear form a(v, u) and linear form L(v) for
# a mixed formulation of Poisson's equation with BDM
# (Brezzi-Douglas-Marini) elements.
#
from ufl import (Coefficient, FiniteElement, TestFunctions, TrialFunctions,
                 div, dot, dx, triangle)

cell = triangle
BDM1 = FiniteElement("Brezzi-Douglas-Marini", cell, 1)
DG0 = FiniteElement("Discontinuous Lagrange", cell, 0)

element = BDM1 * DG0

(tau, w) = TestFunctions(element)
(sigma, u) = TrialFunctions(element)

f = Coefficient(DG0)

a = (dot(tau, sigma) - div(tau) * u + w * div(sigma)) * dx
L = w * f * dx
