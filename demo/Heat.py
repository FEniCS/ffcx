# Copyright (C) 2005-2009 Anders Logg
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
# The bilinear form a(v, u1) and linear form L(v) for
# one backward Euler step with the heat equation.
#
from ufl import (Coefficient, Constant, FiniteElement, TestFunction,
                 TrialFunction, dot, dx, grad, triangle)

cell = triangle
element = FiniteElement("Lagrange", cell, 1)

v = TestFunction(element)  # Test function
u1 = TrialFunction(element)  # Value at t_n
u0 = Coefficient(element)      # Value at t_n-1
c = Coefficient(element)      # Heat conductivity
f = Coefficient(element)      # Heat source
k = Constant(cell)         # Time step

a = v * u1 * dx + k * c * dot(grad(v), grad(u1)) * dx
L = v * u0 * dx + k * v * f * dx
