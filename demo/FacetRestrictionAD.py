# Copyright (C) 2010 Garth N. Wells
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
"""Facet restriction demo."""

import basix.ufl
from ufl import (
    Coefficient,
    FunctionSpace,
    Mesh,
    TestFunction,
    TrialFunction,
    avg,
    derivative,
    dS,
    dx,
    grad,
    inner,
)

element = basix.ufl.element("Discontinuous Lagrange", "triangle", 1)
domain = Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
space = FunctionSpace(domain, element)

v = TestFunction(space)
w = Coefficient(space)
L = inner(grad(w), grad(v)) * dx - inner(avg(grad(w)), avg(grad(v))) * dS

u = TrialFunction(space)
a = derivative(L, w, u)
