# Copyright (C) 2009 Peter Brune
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
"""Normals demo.

This example demonstrates how to use the facet normals
Merely project the normal onto a vector section.
"""

import basix.ufl
from ufl import FacetNormal, FunctionSpace, Mesh, TestFunction, TrialFunction, ds, inner, triangle

cell = triangle

element = basix.ufl.element("Lagrange", cell.cellname(), 1, shape=(2,))
domain = Mesh(basix.ufl.element("Lagrange", cell.cellname(), 1, shape=(2,)))
space = FunctionSpace(domain, element)

n = FacetNormal(domain)

v = TrialFunction(space)
u = TestFunction(space)

a = inner(v, u) * ds
L = inner(n, u) * ds
