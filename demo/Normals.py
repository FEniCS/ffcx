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
#
# This example demonstrates how to use the facet normals
# Merely project the normal onto a vector section.
import basix.ufl
from ufl import FacetNormal, TestFunction, TrialFunction, dot, ds, triangle

cell = triangle

element = basix.ufl.element("Lagrange", cell.cellname(), 1, shape=(2, ))

n = FacetNormal(cell)

v = TrialFunction(element)
u = TestFunction(element)

a = dot(v, u) * ds
L = dot(n, u) * ds
