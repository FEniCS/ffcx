# Copyright (C) 2006-2007 Anders Logg
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
# The bilinear form a(v, u) and linear form L(v) for
# Poisson's equation with Neumann boundary conditions.
from ufl import (Coefficient, TestFunction, TrialFunction, VectorElement, ds,
                 dx, grad, inner, triangle)

element = VectorElement("Lagrange", triangle, 1)

v = TestFunction(element)
u = TrialFunction(element)
f = Coefficient(element)
g = Coefficient(element)

a = inner(grad(v), grad(u)) * dx
L = inner(v, f) * dx + inner(v, g) * ds
