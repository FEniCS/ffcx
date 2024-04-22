# Copyright (C) 2023 Igor A. Baratta
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Mass action demo."""

import basix
import ufl

P = 3
cell_type = basix.CellType.hexahedron
# create element with tensor product order
element = basix.ufl.wrap_element(
    basix.create_tp_element(basix.ElementFamily.P, cell_type, P, basix.LagrangeVariant.gll_warped)
)

coords = basix.ufl.element(basix.ElementFamily.P, cell_type, 1, shape=(3,))
mesh = ufl.Mesh(coords)
V = ufl.FunctionSpace(mesh, element)
x = ufl.SpatialCoordinate(mesh)

v = ufl.TestFunction(V)
u = ufl.TrialFunction(V)
a = ufl.inner(u, v) * ufl.dx

w = ufl.Coefficient(V)
L = ufl.action(a, w)
