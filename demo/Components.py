# Copyright (C) 2011 Garth N. Wells
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
"""Components demo.

This example demonstrates how to create vectors component-wise.
"""

import basix.ufl
from ufl import Coefficient, FunctionSpace, Mesh, TestFunction, as_vector, dx, inner

element = basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,))
domain = Mesh(element)
space = FunctionSpace(domain, element)

v = TestFunction(space)
f = Coefficient(space)

# Create vector
v0 = as_vector([v[0], v[1], 0.0])

# Use created vector in linear form
L = inner(f, v0) * dx
