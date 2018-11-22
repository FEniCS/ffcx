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


def fiat_coordinate_mapping(L, cellname, gdim, ref_coord_symbol="Y"):

    # Code used in evaluatebasis[|derivatives]
    x = L.Symbol("x")
    Y = L.Symbol(ref_coord_symbol)
    coord_dofs = L.Symbol("coordinate_dofs")

    if cellname == "interval":
        J = L.Symbol("J")
        detJ = L.Symbol("detJ")
        if gdim == 1:
            code = [
                L.Comment("Get coordinates and map to the reference (FIAT) element"),
                L.ArrayDecl("double", Y, 1,
                            [(2 * x[0] - coord_dofs[0] - coord_dofs[1]) / J[0]])
            ]
        elif gdim == 2:
            code = [
                L.Comment("Get coordinates and map to the reference (FIAT) element"),
                L.ArrayDecl("double", Y, 1, [
                    2 * (L.Sqrt(
                        L.Call("pow", (x[0] - coord_dofs[0], 2))
                        + L.Call("pow", (x[1] - coord_dofs[1], 2)))) / detJ - 1.0
                ])
            ]
        elif gdim == 3:
            code = [
                L.Comment("Get coordinates and map to the reference (FIAT) element"),
                L.ArrayDecl("double", Y, 1, [
                    2 * (L.Sqrt(
                        L.Call("pow", (x[0] - coord_dofs[0], 2))
                        + L.Call("pow", (x[1] - coord_dofs[1], 2))
                        + L.Call("pow", (x[2] - coord_dofs[2], 2)))) / detJ - 1.0
                ])
            ]
        else:
            raise RuntimeError("Cannot compute interval with gdim: %d" % gdim)
    elif cellname == "triangle":
        if gdim == 2:
            C0 = L.Symbol("C0")
            C1 = L.Symbol("C1")
            J = L.Symbol("J")
            detJ = L.Symbol("detJ")
            code = [
                L.Comment("Compute constants"),
                L.VariableDecl("const double", C0, coord_dofs[2] + coord_dofs[4]),
                L.VariableDecl("const double", C1, coord_dofs[3] + coord_dofs[5]),
                L.Comment("Get coordinates and map to the reference (FIAT) element"),
                L.ArrayDecl("double", Y, 2,
                            [(J[1] * (C1 - 2.0 * x[1]) + J[3] * (2.0 * x[0] - C0)) / detJ,
                             (J[0] * (2.0 * x[1] - C1) + J[2] * (C0 - 2.0 * x[0])) / detJ])
            ]
        elif gdim == 3:
            K = L.Symbol("K")
            code = [
                L.Comment("P_FFC = J^dag (p - b), P_FIAT = 2*P_FFC - (1, 1)"),
                L.ArrayDecl("double", Y, 2, [
                    2 * (K[0] * (x[0] - coord_dofs[0]) + K[1] * (x[1] - coord_dofs[1])
                         + K[2] * (x[2] - coord_dofs[2])) - 1.0,
                    2 * (K[3] * (x[0] - coord_dofs[0]) + K[4] * (x[1] - coord_dofs[1])
                         + K[5] * (x[2] - coord_dofs[2])) - 1.0
                ])
            ]
        else:
            raise RuntimeError("Cannot compute triangle with gdim: %d" % gdim)
    elif cellname == 'tetrahedron' and gdim == 3:
        C0 = L.Symbol("C0")
        C1 = L.Symbol("C1")
        C2 = L.Symbol("C2")
        J = L.Symbol("J")
        detJ = L.Symbol("detJ")
        d = L.Symbol("d")

        code = [
            L.Comment("Compute constants"),
            L.VariableDecl(
                "const double", C0,
                coord_dofs[9] + coord_dofs[6] + coord_dofs[3] - coord_dofs[0]),
            L.VariableDecl(
                "const double", C1,
                coord_dofs[10] + coord_dofs[7] + coord_dofs[4] - coord_dofs[1]),
            L.VariableDecl(
                "const double", C2,
                coord_dofs[11] + coord_dofs[8] + coord_dofs[5] - coord_dofs[2]),
            L.Comment("Compute subdeterminants"),
            L.ArrayDecl("const double", d, 9, [
                J[4] * J[8] - J[5] * J[7], J[5] * J[6] - J[3] * J[8], J[3] * J[7] - J[4] * J[6],
                J[2] * J[7] - J[1] * J[8], J[0] * J[8] - J[2] * J[6], J[1] * J[6] - J[0] * J[7],
                J[1] * J[5] - J[2] * J[4], J[2] * J[3] - J[0] * J[5], J[0] * J[4] - J[1] * J[3]
            ]),
            L.Comment("Get coordinates and map to the reference (FIAT) element"),
            L.ArrayDecl("double", Y, 3, [(d[0] * (2.0 * x[0] - C0) + d[3]
                                          * (2.0 * x[1] - C1) + d[6] * (2.0 * x[2] - C2)) / detJ,
                                         (d[1] * (2.0 * x[0] - C0) + d[4]
                                          * (2.0 * x[1] - C1) + d[7] * (2.0 * x[2] - C2)) / detJ,
                                         (d[2] * (2.0 * x[0] - C0) + d[5]
                                          * (2.0 * x[1] - C1) + d[8] * (2.0 * x[2] - C2)) / detJ])
        ]
    else:
        raise RuntimeError("Cannot compute %s with gdim: %d" % (cellname, gdim))

    return code
