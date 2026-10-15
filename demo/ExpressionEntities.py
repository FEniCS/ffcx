# Copyright (C) 2026 Jørgen S. Dokken
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
# along with FFC. If not, see <http://www.gnu.org/licenses/>.
"""Expressions on mesh entities demo.

Defines Expressions evaluated on cells, facets, ridges and peaks. The
dimension of the evaluation points selects the entity: points on a
vertex have zero columns, which the generated code must handle without
a zero-size array.
"""

import basix.ufl
import numpy as np
from ufl import Coefficient, FunctionSpace, Mesh, SpatialCoordinate, grad

# A triangle mesh, whose ridges are vertices
mesh = Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
V = FunctionSpace(mesh, basix.ufl.element("Lagrange", "triangle", 2))
u = Coefficient(V)
x = SpatialCoordinate(mesh)

cell_points = np.array([[1.0 / 3.0, 1.0 / 3.0]])
facet_points = np.array([[0.25], [0.75]])
vertex_points = np.zeros((1, 0))

# An interval mesh, whose facets are vertices
mesh_1d = Mesh(basix.ufl.element("Lagrange", "interval", 1, shape=(1,)))
V_1d = FunctionSpace(mesh_1d, basix.ufl.element("Lagrange", "interval", 2))
u_1d = Coefficient(V_1d)

# Data on a mesh of the triangle's vertices
point_mesh = Mesh(basix.ufl.element("Lagrange", "point", 0, shape=(2,), discontinuous=True))
W = FunctionSpace(point_mesh, basix.ufl.element("Lagrange", "point", 0))
w = Coefficient(W)

# A tetrahedron mesh, whose peaks are vertices, with data on a mesh of them
mesh_3d = Mesh(basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,)))
V_3d = FunctionSpace(mesh_3d, basix.ufl.element("Lagrange", "tetrahedron", 2))
u_3d = Coefficient(V_3d)
point_mesh_3d = Mesh(basix.ufl.element("Lagrange", "point", 0, shape=(3,), discontinuous=True))
w_3d = Coefficient(FunctionSpace(point_mesh_3d, basix.ufl.element("Lagrange", "point", 0)))

expressions = [
    (u * x[0], cell_points),
    (grad(u), facet_points),
    (grad(u), vertex_points),
    (u_1d.dx(0), vertex_points),
    (w * u * SpatialCoordinate(point_mesh)[1], vertex_points),
    (grad(u_3d), vertex_points),
    (w_3d * u_3d * SpatialCoordinate(point_mesh_3d)[2], vertex_points),
]
