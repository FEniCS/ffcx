# Copyright (C) 2010 Kristian B. Oelgaard
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
"""Spatial coordinates demo.

The bilinear form a(u, v) and linear form L(v) for Poisson's equation where
spatial coordinates are used to define the source and boundary flux terms.
"""

import basix.ufl
from ufl import (
    FunctionSpace,
    Mesh,
    SpatialCoordinate,
    TestFunction,
    TrialFunction,
    ds,
    dx,
    exp,
    grad,
    inner,
    sin,
)

element = basix.ufl.element("Lagrange", "triangle", 2)
domain = Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
space = FunctionSpace(domain, element)

u = TrialFunction(space)
v = TestFunction(space)

x = SpatialCoordinate(domain)
d_x = x[0] - 0.5
d_y = x[1] - 0.5
f = 10.0 * exp(-(d_x * d_x + d_y * d_y) / 0.02)
g = sin(5.0 * x[0])

a = inner(grad(u), grad(v)) * dx
L = inner(f, v) * dx + inner(g, v) * ds
