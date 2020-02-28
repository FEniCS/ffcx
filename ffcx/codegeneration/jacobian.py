# Copyright (C) 2009-2017 Anders Logg and Martin Sandve Alnæs, Chris Richardson
#
# This file is part of FFCX.(https://www.fenicsproject.org)
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
