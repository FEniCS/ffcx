# -*- coding: utf-8 -*-
# Copyright (C) 2009-2017 Anders Logg and Martin Sandve Alnæs, Chris Richardson
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

# Note: Much of the code in this file is a direct translation
# from the old implementation in FFC, although some improvements
# have been made to the generated code.


def jacobian(L, gdim, tdim, element_cellname):
    J = L.Symbol("J")
    coord_dofs = L.Symbol("coordinate_dofs")
    code = [
        L.Comment("Compute Jacobian"),
        L.ArrayDecl("double", J, (gdim * tdim, )),
        L.Call("compute_jacobian_" + element_cellname + "_" + str(gdim) + "d", (J, coord_dofs))
    ]
    return code


def inverse_jacobian(L, gdim, tdim, element_cellname):
    K = L.Symbol("K")
    J = L.Symbol("J")
    detJ = L.Symbol("detJ")
    code = [
        L.Comment("Compute Inverse Jacobian and determinant"),
        L.ArrayDecl("double", K, (gdim * tdim, )),
        L.VariableDecl("double", detJ),
        L.Call("compute_jacobian_inverse_" + element_cellname + "_" + str(gdim) + "d",
               (K, L.AddressOf(detJ), J))
    ]
    return code


def orientation(L):
    detJ = L.Symbol("detJ")
    cell_orientation = L.Symbol("cell_orientation")
    code = [
        L.Comment("Check orientation and return error code"),
        L.If(L.EQ(cell_orientation, -1), [L.Return(-1)]),
        L.Comment("(If cell_orientation == 1 = down, multiply det(J) by -1)"),
        L.ElseIf(L.EQ(cell_orientation, 1), [L.AssignMul(detJ, -1)])
    ]
    return code
