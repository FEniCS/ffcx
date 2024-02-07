# Copyright (C) 2004-2007 Anders Logg
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
"""1D Poisson demo.

The bilinear form a(u, v) and linear form L(v) for Poisson's equation.
"""

import basix.ufl
from ufl import Coefficient, FunctionSpace, Mesh, TestFunction, TrialFunction, dx, grad, inner

element = basix.ufl.element("Lagrange", "interval", 1)
domain = Mesh(basix.ufl.element("Lagrange", "interval", 1, shape=(1,)))
space = FunctionSpace(domain, element)

u = TrialFunction(space)
v = TestFunction(space)
f = Coefficient(space)

a = inner(grad(u), grad(v)) * dx
L = inner(f, v) * dx
