# Copyright (C) 2023 Igor A. Baratta
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import numpy as np
import basix
import ufl

P = 3
cell_type = basix.CellType.hexahedron
dof_ordering = np.array([0, 16, 4, 20, 1, 17, 5, 21, 32, 48, 8, 12, 2,
                         3, 24, 28, 18, 19, 36, 52, 6, 7, 22, 23, 33, 49,
                         9, 13, 25, 29, 37, 53, 40, 56, 44, 60, 34, 50, 35,
                         51, 10, 14, 11, 15, 26, 30, 27, 31, 38, 54, 39, 55,
                         41, 57, 45, 61, 42, 58, 46, 62, 43, 59, 47, 63])

# create element with tensor product order
element_tp = basix.create_element(basix.ElementFamily.P, cell_type,
                                  P, basix.LagrangeVariant.gll_warped,
                                  dof_ordering=dof_ordering)
element = basix.ufl._BasixElement(element_tp)

coords = basix.ufl.element(basix.ElementFamily.P, cell_type, 1, shape=(3, ))
mesh = ufl.Mesh(coords)
V = ufl.FunctionSpace(mesh, element)
x = ufl.SpatialCoordinate(mesh)

v = ufl.TestFunction(V)
u = ufl.TrialFunction(V)
a = x[0] * ufl.inner(u, v) * ufl.dx

w = ufl.Coefficient(V)
L = ufl.action(a, w)
