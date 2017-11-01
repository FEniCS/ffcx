# -*- coding: utf-8 -*-
# Copyright (C) 2009-2017 Anders Logg and Martin Sandve Alnæs, Chris Richardson
#
# This file is part of UFLACS.
#
# UFLACS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# UFLACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with UFLACS. If not, see <http://www.gnu.org/licenses/>.


# Note: Much of the code in this file is a direct translation
# from the old implementation in FFC, although some improvements
# have been made to the generated code.

def jacobian(L, gdim, tdim, element_cellname):
    J = L.Symbol("J")
    coordinate_dofs = L.Symbol("coordinate_dofs")
    code = [L.Comment("Compute Jacobian"),
            L.ArrayDecl("double", J, (gdim*tdim,)),
            L.Call("compute_jacobian_"+element_cellname+"_"+str(gdim)+"d",
                   (J, coordinate_dofs))]
    return code

def inverse_jacobian(L, gdim, tdim, element_cellname):
    K = L.Symbol("K")
    J = L.Symbol("J")
    detJ = L.Symbol("detJ")
    code = [L.Comment("Compute Inverse Jacobian and determinant"),
            L.ArrayDecl("double", K, (gdim*tdim,)),
            L.VariableDecl("double", detJ),
            L.Call("compute_jacobian_inverse_"+element_cellname+"_"+str(gdim)+"d",(K, detJ, J))]
    return code

def orientation(L):
    detJ = L.Symbol("detJ")
    cell_orientation = L.Symbol("cell_orientation")
    code = [L.Comment("Check orientation"),
            L.If(L.EQ(cell_orientation, -1),
                 [L.Throw("std::runtime_error", "cell orientation must be defined (not -1)")]),
            L.Comment("(If cell_orientation == 1 = down, multiply det(J) by -1)"),
            L.ElseIf(L.EQ(cell_orientation, 1),
                     [L.AssignMul(detJ , -1)])]
    return code

def fiat_coordinate_mapping(L, cellname, gdim, ref_coord_symbol="Y"):

    # Code used in evaluatebasis[|derivatives]
    x = L.Symbol("x")
    Y = L.Symbol(ref_coord_symbol)
    coordinate_dofs = L.Symbol("coordinate_dofs")

    if cellname == "interval":
        J = L.Symbol("J")
        detJ = L.Symbol("detJ")
        if gdim == 1:
            code = [L.Comment("Get coordinates and map to the reference (FIAT) element"),
                    L.ArrayDecl("double", Y, 1, [(2*x[0] - coordinate_dofs[0] - coordinate_dofs[1])/J[0]])]
        elif gdim == 2:
            code = [L.Comment("Get coordinates and map to the reference (FIAT) element"),
                    L.ArrayDecl("double", Y, 1, [2*(L.Sqrt(L.Call("std::pow", (x[0] - coordinate_dofs[0], 2)) + L.Call("std::pow", (x[1] - coordinate_dofs[1], 2))))/detJ - 1.0])]
        elif gdim == 3:
            code = [L.Comment("Get coordinates and map to the reference (FIAT) element"),
                    L.ArrayDecl("double", Y, 1, [2*(L.Sqrt(L.Call("std::pow", (x[0] - coordinate_dofs[0], 2)) + L.Call("std::pow", (x[1] - coordinate_dofs[1], 2)) + L.Call("std::pow", (x[2] - coordinate_dofs[2], 2))))/ detJ - 1.0])]
        else:
            error("Cannot compute interval with gdim: %d" % gdim)
    elif cellname == "triangle":
        if gdim == 2:
            C0 = L.Symbol("C0")
            C1 = L.Symbol("C1")
            J = L.Symbol("J")
            detJ = L.Symbol("detJ")
            code = [L.Comment("Compute constants"),
                    L.VariableDecl("const double", C0, coordinate_dofs[2] + coordinate_dofs[4]),
                    L.VariableDecl("const double", C1, coordinate_dofs[3] + coordinate_dofs[5]),
                    L.Comment("Get coordinates and map to the reference (FIAT) element"),
                    L.ArrayDecl("double", Y, 2, [(J[1]*(C1 - 2.0*x[1]) + J[3]*(2.0*x[0] - C0)) / detJ,
                                                 (J[0]*(2.0*x[1] - C1) + J[2]*(C0 - 2.0*x[0])) / detJ])]
        elif gdim == 3:
            K = L.Symbol("K")
            code = [L.Comment("P_FFC = J^dag (p - b), P_FIAT = 2*P_FFC - (1, 1)"),
                    L.ArrayDecl("double", Y, 2, [2*(K[0]*(x[0] - coordinate_dofs[0])
                                                    + K[1]*(x[1] - coordinate_dofs[1])
                                                    + K[2]*(x[2] - coordinate_dofs[2])) - 1.0,
                                                 2*(K[3]*(x[0] - coordinate_dofs[0])
                                                    + K[4]*(x[1] - coordinate_dofs[1])
                                                    + K[5]*(x[2] - coordinate_dofs[2])) - 1.0])]
        else:
            error("Cannot compute triangle with gdim: %d" % gdim)
    elif cellname == 'tetrahedron' and gdim == 3:
        C0 = L.Symbol("C0")
        C1 = L.Symbol("C1")
        C2 = L.Symbol("C2")
        J = L.Symbol("J")
        detJ = L.Symbol("detJ")
        d = L.Symbol("d")

        code = [L.Comment("Compute constants"),
                L.VariableDecl("const double", C0, coordinate_dofs[9]  + coordinate_dofs[6] + coordinate_dofs[3] - coordinate_dofs[0]),
                L.VariableDecl("const double", C1, coordinate_dofs[10] + coordinate_dofs[7] + coordinate_dofs[4] - coordinate_dofs[1]),
                L.VariableDecl("const double", C2, coordinate_dofs[11] + coordinate_dofs[8] + coordinate_dofs[5] - coordinate_dofs[2]),
                L.Comment("Compute subdeterminants"),
                L.ArrayDecl("const double", d, 9, [J[4]*J[8] - J[5]*J[7],
                                                   J[5]*J[6] - J[3]*J[8],
                                                   J[3]*J[7] - J[4]*J[6],
                                                   J[2]*J[7] - J[1]*J[8],
                                                   J[0]*J[8] - J[2]*J[6],
                                                   J[1]*J[6] - J[0]*J[7],
                                                   J[1]*J[5] - J[2]*J[4],
                                                   J[2]*J[3] - J[0]*J[5],
                                                   J[0]*J[4] - J[1]*J[3]]),
                L.Comment("Get coordinates and map to the reference (FIAT) element"),
                L.ArrayDecl("double", Y, 3, [(d[0]*(2.0*x[0] - C0) + d[3]*(2.0*x[1] - C1) + d[6]*(2.0*x[2] - C2)) / detJ,
                                             (d[1]*(2.0*x[0] - C0) + d[4]*(2.0*x[1] - C1) + d[7]*(2.0*x[2] - C2)) / detJ,
                                             (d[2]*(2.0*x[0] - C0) + d[5]*(2.0*x[1] - C1) + d[8]*(2.0*x[2] - C2)) / detJ])]
    else:
        error("Cannot compute %s with gdim: %d" % (cellname, gdim))

    return code


def _mapping_transform(L, data, dof_data, values, offset, width=1):

    code = []
    tdim = data["topological_dimension"]
    gdim = data["geometric_dimension"]
    num_components = dof_data["num_components"]
    mapping = dof_data["mapping"]

    # Apply transformation if applicable.
    if mapping == "affine":
        return code

    # Get temporary values before mapping.
    tmp_ref = L.Symbol("tmp_ref")
    code += [L.ArrayDecl("const double", tmp_ref, num_components,
                [values[i*width + offset] for i in range(num_components)])]

    # Define symbols for Jacobian, inverse and determinant
    J = L.Symbol("J")
    J = L.FlattenedArray(J, dims=(gdim, tdim))
    detJ = L.Symbol("detJ")
    K = L.Symbol("K")
    K = L.FlattenedArray(K, dims=(tdim, gdim))

    if mapping == "contravariant piola":
        code += [L.Comment("Using contravariant Piola transform to map values back to the physical element")]
        assert num_components == tdim
        for i in range(gdim):
            inner = sum(J[i, j]*tmp_ref[j] for j in range(tdim))
            code += [L.Assign(values[i*width + offset], inner/detJ)]

    elif mapping == "covariant piola":
        code += [L.Comment("Using covariant Piola transform to map values back to the physical element")]
        assert num_components == tdim
        for i in range(gdim):
            inner = sum(K[j, i]*tmp_ref[j] for j in range(tdim))
            code += [L.Assign(values[i*width + offset], inner)]

    elif mapping == "double covariant piola":
        code += [L.Comment("Using double covariant Piola transform to map values back to the physical element")]
        assert num_components == tdim**2
        for p in range(num_components):
            # unflatten the indices
            i = p // tdim
            l = p % tdim  # noqa: E741
            # g_il = K_ji G_jk K_kl
            inner = sum(K[j, i]*tmp_ref[j*tdim + k]*K[k, l]
                        for j in range(tdim) for k in range(tdim))

            code += [L.Assign(values[p*width + offset], inner)]

    elif mapping == "double contravariant piola":
        code += [L.Comment("Using double contravariant Piola transform to map values back to the physical element")]
        assert num_components == tdim**2
        for p in range(num_components):
            # unflatten the indices
            i = p // tdim
            l = p % tdim  # noqa: E741
            inner = sum(J[i, j]*tmp_ref[j*tdim + k]*J[l, k]
                        for j in range(tdim) for k in range(tdim))

            code += [L.Assign(values[p*width + offset], inner/(detJ*detJ))]
    else:
        error("Unknown mapping: %s" % mapping)
    return code
