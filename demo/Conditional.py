# Copyright (C) 2010-2011 Kristian B. Oelgaard
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
"""Conditional demo.

Illustration on how to use Conditional to define a source term.
"""

import basix.ufl
from ufl import (
    And,
    Constant,
    FunctionSpace,
    Mesh,
    Not,
    Or,
    SpatialCoordinate,
    TestFunction,
    conditional,
    dx,
    ge,
    gt,
    inner,
    le,
    lt,
)

element = basix.ufl.element("Lagrange", "triangle", 2)
domain = Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
space = FunctionSpace(domain, element)

v = TestFunction(space)
g = Constant(domain)

x = SpatialCoordinate(domain)
c0 = conditional(le((x[0] - 0.33) ** 2 + (x[1] - 0.67) ** 2, 0.015), -1.0, 5.0)
c = conditional(le((x[0] - 0.33) ** 2 + (x[1] - 0.67) ** 2, 0.025), c0, 0.0)

t0 = And(ge(x[0], 0.55), le(x[0], 0.95))
t1 = Or(lt(x[1], 0.05), gt(x[1], 0.45))
t2 = And(t0, Not(t1))
t = conditional(And(ge(x[1] - x[0] - 0.05 + 0.55, 0.0), t2), -1.0, 0.0)

k = conditional(gt(1, 0), g, g + 1)

f = c + t + k

L = inner(f, v) * dx
