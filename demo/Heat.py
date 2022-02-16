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
# The bilinear form a(v, u1) and linear form L(v) for
# one backward Euler step with the heat equation.
#
# Compile this form with FFCx: ffcx Heat.ufl
from ufl import (Coefficient, Constant, FiniteElement, TestFunction,
                 TrialFunction, dx, grad, inner, triangle)

element = FiniteElement("Lagrange", triangle, 1)

u1 = TrialFunction(element)    # Value at t_n
u0 = Coefficient(element)      # Value at t_n-1
v = TestFunction(element)     # Test function
c = Coefficient(element)      # Heat conductivity
f = Coefficient(element)      # Heat source
k = Constant(triangle)        # Time step

a = u1 * v * dx + k * c * inner(grad(u1), grad(v)) * dx
L = u0 * v * dx + k * f * v * dx
