# Copyright (C) 2009 Kristian B. Oelgaard
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
"""Metadata demo.

Test form for metadata.
"""

import basix.ufl
from ufl import (
    Coefficient,
    Constant,
    FunctionSpace,
    Mesh,
    TestFunction,
    TrialFunction,
    dx,
    grad,
    inner,
)

element = basix.ufl.element("Lagrange", "triangle", 1)
vector_element = basix.ufl.element("Lagrange", "triangle", 1, shape=(2,))
domain = Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
space = FunctionSpace(domain, element)
vector_space = FunctionSpace(domain, vector_element)

u = TrialFunction(space)
v = TestFunction(space)
c = Coefficient(vector_space)
c2 = Constant(domain)

# Terms on the same subdomain using different quadrature degree
a = (
    inner(grad(u), grad(v)) * dx(0, degree=8)
    + inner(c, c) * inner(grad(u), grad(v)) * dx(1, degree=4)
    + inner(c, c) * inner(grad(u), grad(v)) * dx(1, degree=2)
    + inner(grad(u), grad(v)) * dx(1, degree=-1)
)

L = inner(c2, v) * dx(0, metadata={"precision": 1})
