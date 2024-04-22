# Copyright (c) 2005-2007 Anders Logg
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
"""Stabilised Stokes demo.

The bilinear form a(u, v) and Linear form L(v) for the Stokes
equations using a mixed formulation (equal-order stabilized).
"""

import basix.ufl
from ufl import (
    Coefficient,
    FunctionSpace,
    Mesh,
    TestFunctions,
    TrialFunctions,
    div,
    dot,
    dx,
    grad,
    inner,
)

vector = basix.ufl.element("Lagrange", "triangle", 1, shape=(2,))
scalar = basix.ufl.element("Lagrange", "triangle", 1)
system = basix.ufl.mixed_element([vector, scalar])
domain = Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
system_space = FunctionSpace(domain, system)
scalar_space = FunctionSpace(domain, scalar)
vector_space = FunctionSpace(domain, vector)

(u, p) = TrialFunctions(system_space)
(v, q) = TestFunctions(system_space)

f = Coefficient(vector_space)
h = Coefficient(scalar_space)

beta = 0.2
delta = beta * h * h

a = (inner(grad(u), grad(v)) - div(v) * p + div(u) * q + delta * dot(grad(p), grad(q))) * dx
L = dot(f, v + delta * grad(q)) * dx
