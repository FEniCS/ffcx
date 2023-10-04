# Copyright (C) 2023 Chris Richardson
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
# The bilinear form a(u, v) and linear form L(v) for
# Poisson's equation using bilinear elements on bilinear mesh geometry.
import basix.ufl
from ufl import (Coefficient, FunctionSpace, Mesh, TestFunction, TrialFunction,
                 dx, grad, inner)

coords = basix.ufl.element("P", "triangle", 2, shape=(2, ))
mesh = Mesh(coords)
dx = dx(mesh)

element = basix.ufl.element("P", mesh.ufl_cell().cellname(), 2)
space = FunctionSpace(mesh, element)

u = TrialFunction(space)
v = TestFunction(space)
f = Coefficient(space)

# Test literal complex number in form
k = 3.213 + 1.023j

a = k * inner(grad(u), grad(v)) * dx
L = inner(k * f, v) * dx
