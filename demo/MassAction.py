# Copyright (C) 2023 Igor A. Baratta
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import ufl
import basix
import numpy as np

P = 3
cell_type = basix.CellType.hexahedron
element = basix.create_element(basix.ElementFamily.P, cell_type,
                               P, basix.LagrangeVariant.gll_warped)

# create element with tensor product order
factors = element.get_tensor_product_representation()[0]
perm = factors[1]

# create element with tensor product order
element_tp = basix.create_element(basix.ElementFamily.P, cell_type,
                                  P, basix.LagrangeVariant.gll_warped,
                                  dof_ordering=np.argsort(perm))
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
