# Copyright (C) 2007 Anders Logg
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
# Test form for scalar and vector constants.
from ufl import (Coefficient, Constant, FiniteElement, TestFunction,
                 TrialFunction, VectorConstant, dot, dx, grad, inner, triangle)

cell = triangle
element = FiniteElement("Lagrange", cell, 1)

v = TestFunction(element)
u = TrialFunction(element)
f = Coefficient(element)

c = Constant(cell)
d = VectorConstant(cell)

a = c * dot(grad(v), grad(u)) * dx
L = inner(d, grad(v)) * dx
