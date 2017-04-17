# -*- coding: utf-8 -*-
# Copyright (C) 2009-2017 Anders Logg and Martin Sandve Alnæs
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


from collections import defaultdict
import numpy
from six import string_types

from ufl import product
from ffc.uflacs.backends.ufc.generator import ufc_generator
from ffc.uflacs.backends.ufc.utils import generate_return_new_switch, generate_return_int_switch, generate_error

from ffc.uflacs.elementtables import clamp_table_small_numbers
from ffc.uflacs.backends.ufc.evaluatebasis import generate_evaluate_reference_basis
from ffc.uflacs.backends.ufc.evaluatebasis import _generate_compute_basisvalues
from ffc.uflacs.backends.ufc.evalderivs import generate_evaluate_reference_basis_derivatives
from ffc.uflacs.backends.ufc.evalderivs import _generate_combinations

# FIXME: Stop depending on legacy code
from ffc.cpp import indent
# from ffc.evaluatebasis import _evaluate_basis
# from ffc.evaluatebasis import _evaluate_basis_all
from ffc.evaluatebasisderivatives import _evaluate_basis_derivatives
from ffc.evaluatebasisderivatives import _evaluate_basis_derivatives_all
# from ffc.interpolatevertexvalues import interpolate_vertex_values
from ffc.evaluatedof import evaluate_dof_and_dofs
# from ffc.evaluatedof import affine_weights
from ufl.permutation import build_component_numbering
from ffc.utils import pick_first

index_type = "std::size_t"

def affine_weights(dim):
    "Compute coefficents for mapping from reference to physical element"

    if dim == 1:
        return lambda x: (1.0 - x[0], x[0])
    elif dim == 2:
        return lambda x: (1.0 - x[0] - x[1], x[0], x[1])
    elif dim == 3:
        return lambda x: (1.0 - x[0] - x[1] - x[2], x[0], x[1], x[2])

def _change_variables(L, mapping, gdim, tdim, offset):
    """Generate code for mapping function values according to
    'mapping' and offset.

    The basics of how to map a field from a physical to the reference
    domain. (For the inverse approach -- see interpolatevertexvalues)

    Let g be a field defined on a physical domain T with physical
    coordinates x. Let T_0 be a reference domain with coordinates
    X. Assume that F: T_0 -> T such that

      x = F(X)

    Let J be the Jacobian of F, i.e J = dx/dX and let K denote the
    inverse of the Jacobian K = J^{-1}. Then we (currently) have the
    following four types of mappings:

    'affine' mapping for g:

      G(X) = g(x)

    For vector fields g:

    'contravariant piola' mapping for g:

      G(X) = det(J) K g(x)   i.e  G_i(X) = det(J) K_ij g_j(x)

    'covariant piola' mapping for g:

      G(X) = J^T g(x)          i.e  G_i(X) = J^T_ij g(x) = J_ji g_j(x)

    'double covariant piola' mapping for g:

      G(X) = J^T g(x) J     i.e. G_il(X) = J_ji g_jk(x) J_kl

    'double contravariant piola' mapping for g:

      G(X) = det(J)^2 K g(x) K^T  i.e. G_il(X)=(detJ)^2 K_ij g_jk K_lk

    """

    # meg: Various mappings must be handled both here and in
    # interpolate_vertex_values. Could this be abstracted out?

    values = L.Symbol("vals")

    if mapping == "affine":
        return [values[offset]]
    elif mapping == "contravariant piola":
        # Map each component from physical to reference using inverse
        # contravariant piola
        detJ = L.Symbol("detJ")
        K = L.Symbol("K")
        K = L.FlattenedArray(K, dims=(tdim, gdim))
        w = []
        for i in range(tdim):
            inner = 0.0
            for j in range(gdim):
                inner += values[j + offset]*K[i, j]
            w.append(inner*detJ)
        return w

    elif mapping == "covariant piola":
        # Map each component from physical to reference using inverse
        # covariant piola
        detJ = L.Symbol("detJ")
        J = L.Symbol("J")
        J = L.FlattenedArray(J, dims=(gdim, tdim))
        w = []
        for i in range(tdim):
            inner = 0.0
            for j in range(gdim):
                inner += values[j + offset]*J[j, i]
            w.append(inner)
        return w

    elif mapping == "double covariant piola":
        # physical to reference pullback as a covariant 2-tensor
        w = []
        J = L.Symbol("J")
        J = L.FlattenedArray(J, dims=(gdim, tdim))
        for i in range(tdim):
            for l in range(tdim):
                inner = 0.0
                for k in range(gdim):
                    for j in range(gdim):
                        inner += J[j, i] * values[j*tdim + k + offset] * J[k, l]
                w.append(inner)
        return w

    elif mapping == "double contravariant piola":
        # physical to reference using double contravariant piola
        w = []
        K = L.Symbol("K")
        K = L.FlattenedArray(K, dims=(tdim, gdim))
        detJ = L.Symbol("detJ")
        for i in range(tdim):
            for l in range(tdim):
                inner = 0.0
                for k in range(gdim):
                    for j in range(gdim):
                        inner += K[i, j] * values[j*tdim + k + offset] * K[l, k]
                w.append(inner*detJ*detJ)
        return w

    else:
        raise Exception("The mapping (%s) is not allowed" % mapping)


def _generate_multiple_points_body(L, i, dof, mapping, gdim, tdim,
                                   offset=0):
    "Generate c++ for-loop for multiple points (integral bodies)"

    result = L.Symbol("result")
    code = [L.Assign(result, 0.0)]
    points = list(dof.keys())
    n = len(points)

    # Get number of tokens per point
    tokens = [dof[x] for x in points]
    len_tokens = pick_first([len(t) for t in tokens])

    # Declare points
    #    points = format["list"]([format["list"](x) for x in points])
    X_i = L.Symbol("X_%d"%i)
    code += [L.ArrayDecl("double", X_i, [n, tdim], points)]

    # Declare components
    components = [[c[0] for (w, c) in token] for token in tokens]
    #    components = format["list"]([format["list"](c) for c in components])
    D_i = L.Symbol("D_%d"%i)
    code += [L.ArrayDecl("int", D_i, [n, len_tokens], components)]

    # Declare weights
    weights = [[w for (w, c) in token] for token in tokens]
    #    weights = format["list"]([format["list"](w) for w in weights])
    W_i = L.Symbol("W_%d"%i)
    code += [L.ArrayDecl("double", W_i, [n, len_tokens], weights)]

    # Declare copy variable:
    copy_i = L.Symbol("copy_%d"%i)
    code += [L.ArrayDecl("double", copy_i, tdim)]

    # Add loop over points
    code += [L.Comment("Loop over points")]

    # Map the points from the reference onto the physical element
    r = L.Symbol("r")
    w0 = L.Symbol("w0")
    w1 = L.Symbol("w1")
    w2 = L.Symbol("w2")
    w3 = L.Symbol("w3")
    d = L.Symbol("d")
    y = L.Symbol("y")
    coordinate_dofs = L.Symbol("coordinate_dofs")
    if tdim == 1:
        lines_r = [L.Comment("Evaluate basis functions for affine mapping"),
                   L.VariableDecl("const double", w0, 1 - X_i*d[r][0]),
                   L.VariableDecl("const double", w1, X_i*d[r][0]),
                   L.Comment("Compute affine mapping y = F(X)")]
        for j in range(gdim):
            lines_r += [L.Assign(y[j],
                                 w0*coordinate_dofs[j] + w1*coordinate_dofs[j + gdim])]
    elif tdim == 2:
        lines_r = [L.Comment("Evaluate basis functions for affine mapping"),
                   L.VariableDecl("const double", w0, 1 - X_i*d[r][0] - X_i*d[r][1]),
                   L.VariableDecl("const double", w1, X_i*d[r][0]),
                   L.VariableDecl("const double", w2, X_i*d[r][1]),
                   L.Comment("Compute affine mapping y = F(X)")]
        for j in range(gdim):
            lines_r += [L.Assign(y[j],
                                 w0*coordinate_dofs[j]
                                 + w1*coordinate_dofs[j + gdim]
                                 + w2*coordinate_dofs[j + 2*gdim])]
    elif tdim == 3:
        lines_r = [L.Comment("Evaluate basis functions for affine mapping"),
                   L.VariableDecl("const double", w0, 1 - X_i*d[r][0] - X_i*d[r][1] - X_i*d[r][2]),
                   L.VariableDecl("const double", w1, X_i*d[r][0]),
                   L.VariableDecl("const double", w2, X_i*d[r][1]),
                   L.VariableDecl("const double", w3, X_i*d[r][2]),
                   L.Comment("Compute affine mapping y = F(X)")]
        for j in range(gdim):
            lines_r += [L.Assign(y[j],
                                 w0*coordinate_dofs[j]
                                 + w1*coordinate_dofs[j + gdim]
                                 + w2*coordinate_dofs[j + 2*gdim]
                                 + w3*coordinate_dofs[j + 3*gdim])]

    # Evaluate function at physical point
    lines_r += [L.Comment("Evaluate function at physical point")]
    y = L.Symbol("y")
    vals = L.Symbol("vals")
    c = L.Symbol("c")
    lines_r += [L.Call("f.evaluate",(y, vals, c))]

    # Map function values to the reference element
    lines_r += [L.Comment("Map function to reference element")]
    F = _change_variables(L, mapping, gdim, tdim, offset)
    lines_r += [L.Assign(copy_i[k], F_k)
                for (k, F_k) in enumerate(F)]

    # Add loop over directional components
    lines_r += [L.Comment("Loop over directions")]

    s = L.Symbol("s")
    lines_r += [L.ForRange(s, 0, len_tokens, index_type=index_type,
                       body=[L.AssignAdd(result, copy_i[D_i[r, s]] * W_i[r, s])])]

    # Generate loop over r and add to code.
    code += [L.ForRange(r, 0, n, index_type=index_type, body=lines_r)]
    #    print(L.StatementList(code))
    return (code, result)

def _generate_body(L, i, dof, mapping, gdim, tdim, offset=0):
    "Generate code for a single dof."

    # EnrichedElement is handled by having [None, ..., None] dual basis
    if not dof:
        msg = "evaluate_dof(s) for enriched element not implemented."
        return generate_error(L, msg, parameters["convert_exceptions_to_warnings"])

    points = list(dof.keys())

    # Generate different code if multiple points. (Otherwise ffc
    # compile time blows up.)
    if len(points) > 1:
        return _generate_multiple_points_body(L, i, dof, mapping, gdim, tdim,
                                              offset)

    # Get weights for mapping reference point to physical
    x = points[0]
    w = affine_weights(tdim)(x)

    # Map point onto physical element: y = F_K(x)
    code = []

    y = L.Symbol("y")
    coordinate_dofs = L.Symbol("coordinate_dofs")
    vals = L.Symbol("vals")
    c = L.Symbol("c")
    for j in range(gdim):
        yy = 0.0
        for k in range(tdim + 1):
            yy += w[k]*coordinate_dofs[k*gdim + j]
        code += [L.Assign(y[j], yy)]

    # Evaluate function at physical point
    code += [L.Call("f.evaluate", (vals, y, c))]

    # Map function values to the reference element
    F = _change_variables(L, mapping, gdim, tdim, offset)

    # Simple affine functions deserve special case:
    if len(F) == 1:
        return (code, dof[x][0][0]*F[0])

    # Flatten multi-indices
    (index_map, _) = build_component_numbering([tdim] * len(dof[x][0][1]), ())
    # Take inner product between components and weights
    value = 0.0
    for (w, k) in dof[x]:
        value += w*F[index_map[k]]

    # Return eval code and value
    return (code, value)

def _x_evaluate_dof_and_dofs(L, ir):
    "Generate code for evaluate_dof and evaluate_dof."

    gdim = ir["geometric_dimension"]
    tdim = ir["topological_dimension"]

    # element_cellname = ir["element_cellname"]
    element_cellname = ['interval', 'triangle', 'tetrahedron'][tdim - 1]

    # Enriched element, no dofs defined
    if not any(ir["dofs"]):
        code = []
    else:
        # Declare variable for storing the result and physical coordinates
        code = [L.Comment("Declare variables for result of evaluation")]
        vals = L.Symbol("vals")
        code += [L.ArrayDecl("double", vals, ir["physical_value_size"])]
        code += [L.Comment("Declare variable for physical coordinates")]
        y = L.Symbol("y")
        code += [L.ArrayDecl("double", y, gdim)]

        # Check whether Jacobians are necessary.
        needs_inverse_jacobian = any(["contravariant piola" in m
                                      for m in ir["mappings"]])
        needs_jacobian = any(["covariant piola" in m for m in ir["mappings"]])

        # Intermediate variable needed for multiple point dofs
        needs_temporary = any(len(dof) > 1 for dof in ir["dofs"])
        if needs_temporary:
            result = L.Symbol("result")
            code += [L.VariableDecl("double", result)]

        if needs_jacobian:
            code += jacobian(L, gdim, tdim, element_cellname)

        if needs_inverse_jacobian:
            code += jacobian(L, gdim, tdim, element_cellname)
            code += inverse_jacobian(L, gdim, tdim, element_cellname)
            code += orientation(L)

    # Extract variables
    mappings = ir["mappings"]
    offsets = ir["physical_offsets"]

    # Generate bodies for each degree of freedom
    cases = []
    for (i, dof) in enumerate(ir["dofs"]):
        c, r = _generate_body(L, i, dof, mappings[i], gdim, tdim, offsets[i])
        c += [L.Return(r)]
        cases.append((i, c))

    code += [L.Switch(L.Symbol("i"), cases)]
    code += [L.Return(0.0)]
    return code

    # Construct dict with eval code as keys to remove duplicate eval code
    #    cases_opt = OrderedDict((case[0], []) for case in cases)
    #for i, (evl, res) in enumerate(cases):
    #cases_opt[evl].append((i, res))
    # Combine each case with assignments for evaluate_dofs
    #dofs_code = reqs
    #for evl, results in six.iteritems(cases_opt):
    #    dofs_code += evl + "\n"
    #    for i, res in results:
    #        dofs_code += format["assign"](component(f_values, i), res) + "\n"
    #dofs_code = dofs_code.rstrip("\n")
    #return (dof_code, dofs_code)

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

def fiat_coordinate_mapping(L, cellname, gdim):

    # Code used in evaluatebasis[|derivatives]
    x = L.Symbol("x")
    Y = L.Symbol("Y")
    coordinate_dofs = L.Symbol("coordinate_dofs")

    if cellname == "interval":
        J = L.Symbol("J")
        detJ = L.Symbol("detJ")
        if gdim == 1:
            code = [L.Comment("Get coordinates and map to the reference (FIAT) element"),
                    L.ArrayDecl("double", Y, 1, [(2*x[0] - coordinate_dofs[0] - coordinate_dofs[1])/J[0]])]
        elif gdim == 2:
            code = [L.Comment("Get coordinates and map to the reference (FIAT) element"),
                    L.ArrayDecl("double", Y, 1, [2*(L.Sqrt(L.Call("std::pow",x[0] - coordinate_dofs[0])) + L.Sqrt(L.Call("std::pow", x[1] - coordinate_dofs[1])))/detJ - 1.0])]
        elif gdim == 3:
            code = [L.Comment("Get coordinates and map to the reference (FIAT) element"),
                    L.ArrayDecl("double", Y, 1, [2*(L.Sqrt(L.Call("std::pow", (x[0] - coordinate_dofs[0], 2)) + L.Call("std::pow", (x[1] - coordinate_dofs[1], 2)) + L.Call("std::pow", (x[2] - coordinate_dofs[2], 2)))/ detJ) - 1.0])]
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

def compute_basis_values(L, data, dof_data):
    basisvalues = L.Symbol("basisvalues")
    Y = L.Symbol("Y")
    element_cellname = data["cellname"]
    embedded_degree = dof_data["embedded_degree"]
    num_members = dof_data["num_expansion_members"]
    return _generate_compute_basisvalues(L, basisvalues, Y, element_cellname, embedded_degree, num_members)

def tabulate_coefficients(L, dof_data):
    """This function tabulates the element coefficients that are
    generated by FIAT at compile time."""

    # Get coefficients from basis functions, computed by FIAT at compile time.
    coefficients = dof_data["coeffs"]

    # Initialise return code.
    code = [L.Comment("Table(s) of coefficients")]

    # Get number of members of the expansion set.
    num_mem = dof_data["num_expansion_members"]

    # Generate tables for each component.
    for i, coeffs in enumerate(coefficients):

        # Variable name for coefficients.
        name = L.Symbol("coefficients%d" % i)

        # Generate array of values.
        code += [L.ArrayDecl("static const double", name, num_mem, coeffs)]

    return code

def compute_values(L, data, dof_data):
    """This function computes the value of the basisfunction as the dot product
    of the coefficients and basisvalues."""

    # Initialise return code.
    code = [L.Comment("Compute value(s)")]

    # Get dof data.
    num_components = dof_data["num_components"]
    reference_offset = dof_data["reference_offset"]
    physical_offset = dof_data["physical_offset"]
    offset = reference_offset  # physical_offset # FIXME: Should be physical offset but that breaks tests

    basisvalues = L.Symbol("basisvalues")
    values = L.Symbol("values")
    r = L.Symbol("r")
    lines = []
    if data["reference_value_size"] != 1:
        # Loop number of components.
        for i in range(num_components):
            coefficients = L.Symbol("coefficients%d" % i)
            lines += [L.AssignAdd(values[i+offset], coefficients[r]*basisvalues[r])]
    else:
        coefficients = L.Symbol("coefficients0")
        lines = [L.AssignAdd(L.Dereference(values), coefficients[r]*basisvalues[r])]

    # Get number of members of the expansion set and generate loop.
    num_mem = dof_data["num_expansion_members"]
    code += [L.ForRange(r, 0, num_mem, index_type=index_type, body=lines)]

    tdim = data["topological_dimension"]
    gdim = data["geometric_dimension"]

    # Apply transformation if applicable.
    mapping = dof_data["mapping"]
    if mapping == "affine":
        pass
    elif mapping == "contravariant piola":
        code += [L.Comment("Using contravariant Piola transform to map values back to the physical element")]

        # Get temporary values before mapping.
        tmp_ref = []
        for i in range(num_components):
            tmp_ref.append(L.Symbol("tmp_ref%d" % i))
        code += [L.VariableDecl("const double", tmp_ref[i], values[i + offset])
                 for i in range(num_components)]


        # Create names for inner product.
        basis_col = [tmp_ref[j] for j in range(tdim)]
        J = L.Symbol("J")
        J = L.FlattenedArray(J, dims=(gdim, tdim))
        detJ = L.Symbol("detJ")
        for i in range(gdim):
            # Create Jacobian.
            jacobian_row = [ J[i, j] for j in range(tdim) ]

            # Create inner product and multiply by inverse of Jacobian.
            inner = 0.0
            for a,b in zip(jacobian_row, basis_col):
                inner += a*b
            value = inner/detJ
            code += [L.Assign(values[i + offset], inner/detJ)]

    elif mapping == "covariant piola":
        code += [L.Comment("Using covariant Piola transform to map values back to the physical element")]
        # Get temporary values before mapping.
        tmp_ref = []
        for i in range(num_components):
            tmp_ref.append(L.Symbol("tmp_ref%d" % i))
        code += [L.VariableDecl("const double", tmp_ref[i], values[i + offset])
                 for i in range(num_components)]

        basis_col = [tmp_ref[j] for j in range(tdim)]
        K = L.Symbol("K")
        K = L.FlattenedArray(K, dims=(tdim, gdim))
        for i in range(gdim):
            # Create inverse of Jacobian.
            inv_jacobian_column = [K[j, i] for j in range(tdim)]

            # Create inner product of basis values and inverse of Jacobian.
            inner = 0.0
            for a, b in zip(inv_jacobian_column, basis_col):
                inner += a*b
            code += [L.Assign(values[i + offset], inner)]

    elif mapping == "double covariant piola":
        code += [L.Comment("Using double covariant Piola transform to map values back to the physical element")]
        # Get temporary values before mapping.
        basis_col = []
        for i in range(num_components):
            basis_col.append(L.Symbol("tmp_ref%d" % i))
        code += [L.VariableDecl("const double", basis_col[i], values[i + offset])
                 for i in range(num_components)]

        # value = f_group(f_inner(
        #     [f_inner([f_trans("JINV", j, i, tdim, gdim, None)
        #               for j in range(tdim)],
        #              [basis_col[j * tdim + k] for j in range(tdim)])
        #      for k in range(tdim)],
        #     [f_trans("JINV", k, l, tdim, gdim, None)
        #      for k in range(tdim)]))

        K = L.Symbol("K")
        for p in range(num_components):
            # unflatten the indices
            i = p // tdim
            l = p % tdim
            # g_il = K_ji G_jk K_kl
            acc_list = []
            for k in range(tdim):
                acc = 0.0
                for j in range(tdim):
                    acc += K[j*gdim +i]*basis_col[j*tdim + k]
                acc_list.append(acc)
            inner = 0.0
            for k in range(tdim):
                inner += acc_list[k]*K[k*gdim + l]

            code += [L.Assign(values[p + offset], inner)]

    elif mapping == "double contravariant piola":
        code += [L.Comment("Using double contravariant Piola transform to map values back to the physical element")]

        # Get temporary values before mapping.
        basis_col = []
        for i in range(num_components):
            basis_col.append(L.Symbol("tmp_ref%d" % i))
        code += [L.VariableDecl("const double", basis_col[i], values[i + offset])
                 for i in range(num_components)]

        J = L.Symbol("J")
        detJ = L.Symbol("detJ")
        for p in range(num_components):
            # unflatten the indices
            i = p // tdim
            l = p % tdim

            # g_il = (det J)^(-2) Jij G_jk Jlk
            #         value = f_group(f_inner(
            #             [f_inner([f_trans("J", i, j, tdim, gdim, None)
            #                       for j in range(tdim)],
            #                      [basis_col[j * tdim + k] for j in range(tdim)])
            #              for k in range(tdim)],
            #             [f_trans("J", l, k, tdim, gdim, None) for k in range(tdim)]))

            acc_list = []
            for k in range(tdim):
                acc = 0.0
                for j in range(tdim):
                    acc += J[i*gdim + j]*basis_col[j*tdim + k]
                acc_list.append(acc)
            inner = 0.0
            for k in range(tdim):
                inner += acc_list[k]*J[l*gdim + k]

            code += [L.Assign(values[p + offset], inner/(detJ*detJ))]
    else:
        error("Unknown mapping: %s" % mapping)

    return code

def generate_element_mapping(mapping, i, num_reference_components, tdim, gdim, J, detJ, K):
    # Select transformation to apply
    if mapping == "affine":
        assert num_reference_components == 1
        num_physical_components = 1
        M_scale = 1
        M_row = [1]  # M_row[0] == 1
    elif mapping == "contravariant piola":
        assert num_reference_components == tdim
        num_physical_components = gdim
        M_scale = 1.0 / detJ
        M_row = [J[i, jj] for jj in range(tdim)]
    elif mapping == "covariant piola":
        assert num_reference_components == tdim
        num_physical_components = gdim
        M_scale = 1.0 / detJ
        M_row = [K[jj, i] for jj in range(tdim)]
    elif mapping == "double covariant piola":
        assert num_reference_components == tdim**2
        num_physical_components = gdim**2
        # g_il = K_ji G_jk K_kl = K_ji K_kl G_jk
        i0 = i // tdim  # i in the line above
        i1 = i % tdim   # l ...
        M_scale = 1.0
        M_row = [K[jj,i0]*K[kk,i1] for jj in range(tdim) for kk in range(tdim)]
    elif mapping == "double contravariant piola":
        assert num_reference_components == tdim**2
        num_physical_components = gdim**2
        # g_il = (det J)^(-2) Jij G_jk Jlk = (det J)^(-2) Jij Jlk G_jk
        i0 = i // tdim  # i in the line above
        i1 = i % tdim   # l ...
        M_scale = 1.0 / (detJ*detJ)
        M_row = [J[i0,jj]*J[i1,kk] for jj in range(tdim) for kk in range(tdim)]
    else:
        error("Unknown mapping: %s" % mapping)
    return M_scale, M_row, num_physical_components


class ufc_finite_element(ufc_generator):
    "Each function maps to a keyword in the template. See documentation of ufc_generator."
    def __init__(self):
        ufc_generator.__init__(self, "finite_element")

    def cell_shape(self, L, cell_shape):
        return L.Return(L.Symbol("ufc::shape::" + cell_shape))

    def topological_dimension(self, L, topological_dimension):
        return L.Return(topological_dimension)

    def geometric_dimension(self, L, geometric_dimension):
        return L.Return(geometric_dimension)

    def space_dimension(self, L, space_dimension):
        return L.Return(space_dimension)

    def value_rank(self, L, value_shape):
        return L.Return(len(value_shape))

    def value_dimension(self, L, value_shape):
        return generate_return_int_switch(L, "i", value_shape, 1)

    def value_size(self, L, value_shape):
        return L.Return(product(value_shape))

    def reference_value_rank(self, L, reference_value_shape):
        return L.Return(len(reference_value_shape))

    def reference_value_dimension(self, L, reference_value_shape):
        return generate_return_int_switch(L, "i", reference_value_shape, 1)

    def reference_value_size(self, L, reference_value_shape):
        return L.Return(product(reference_value_shape))

    def degree(self, L, degree):
        return L.Return(degree)

    def family(self, L, family):
        return L.Return(L.LiteralString(family))

    def num_sub_elements(self, L, num_sub_elements):
        return L.Return(num_sub_elements)

    def create_sub_element(self, L, ir):
        classnames = ir["create_sub_element"]
        return generate_return_new_switch(L, "i", classnames, factory=ir["jit"])

    def evaluate_basis(self, L, ir, parameters):
        #        legacy_code = indent(_evaluate_basis(ir["evaluate_basis"]), 4)
        #        print(legacy_code)

        data = ir["evaluate_basis"]

        # FIXME: does this make sense?
        if not data:
            msg = "evaluate_basis is not defined for this element"
            return generate_error(L, msg, parameters["convert_exceptions_to_warnings"])

        # Get the element cell name and geometric dimension.
        element_cellname = data["cellname"]
        gdim = data["geometric_dimension"]
        tdim = data["topological_dimension"]

        # Generate run time code to evaluate an element basisfunction at an
        # arbitrary point. The value(s) of the basisfunction is/are
        # computed as in FIAT as the dot product of the coefficients (computed at compile time)
        # and basisvalues which are dependent on the coordinate and thus have to be computed at
        # run time.

        # The function should work for all elements supported by FIAT, but it remains
        # untested for tensor valued elements.

        # Get code snippets for Jacobian, Inverse of Jacobian and mapping of
        # coordinates from physical element to the FIAT reference element.

        code = jacobian(L, gdim, tdim, element_cellname)
        code += inverse_jacobian(L, gdim, tdim, element_cellname)
        if data["needs_oriented"]:
            code += orientation(L)

        if any((d["embedded_degree"] > 0) for d in data["dofs_data"]):
            code += fiat_coordinate_mapping(L, element_cellname, gdim)

        reference_value_size = data["reference_value_size"]
        code += [L.Comment("Reset values")]
        dof_values = L.Symbol("values")
        if reference_value_size == 1:
            # Reset values as a pointer.
            code += [L.Assign(L.Dereference(dof_values), 0.0)]
        else:
            code += [L.MemZero(dof_values, reference_value_size)]

        # Create code for all basis values (dofs).
        dof_cases = []
        for f, dof_data in enumerate(data["dofs_data"]):
            dof_code = compute_basis_values(L, data, dof_data)
            dof_code += tabulate_coefficients(L, dof_data)
            dof_code += compute_values(L, data, dof_data)
            dof_cases.append((f, dof_code))

        code += [L.Switch(L.Symbol("i"), dof_cases)]
        #        print(L.StatementList(code))
        return code

    def evaluate_basis_all(self, L, ir, parameters):

        data=ir["evaluate_basis"]
        physical_value_size = data["physical_value_size"]
        space_dimension = data["space_dimension"]

        x = L.Symbol("x")
        coordinate_dofs = L.Symbol("coordinate_dofs")
        cell_orientation = L.Symbol("cell_orientation")
        values = L.Symbol("values")

        # Special case where space dimension is one (constant elements).
        if space_dimension == 1:
            code = [L.Comment("Element is constant, calling evaluate_basis."),
                    L.Call("evaluate_basis",
                           (0, values, x, coordinate_dofs, cell_orientation))]
            return code

        r = L.Symbol("r")
        dof_values = L.Symbol("dof_values")
        if physical_value_size == 1:
            code = [ L.Comment("Helper variable to hold value of a single dof."),
                     L.VariableDecl("double", dof_values, 0.0),
                     L.Comment("Loop dofs and call evaluate_basis"),
                     L.ForRange(r, 0, space_dimension, index_type=index_type,
                                body=[L.Call("evaluate_basis",
                                             (r, L.AddressOf(dof_values), x,
                                              coordinate_dofs, cell_orientation)),
                                      L.Assign(values[r], dof_values)]
                               )
                   ]
        else:
            s = L.Symbol("s")
            code = [L.Comment("Helper variable to hold values of a single dof."),
                    L.ArrayDecl("double", dof_values, physical_value_size, 0.0),
                    L.Comment("Loop dofs and call evaluate_basis"),
                    L.ForRange(r, 0, space_dimension, index_type=index_type,
                               body=[L.Call("evaluate_basis",
                                             (r, dof_values, x,
                                              coordinate_dofs, cell_orientation)),
                                     L.ForRange(s, 0, physical_value_size,
                                                index_type=index_type,
                                                body=[L.Assign(values[r*physical_value_size+s], dof_values[s])])
                                    ]
                              )
                   ]

        return code

    def evaluate_basis_derivatives(self, L, ir, parameters):
        # FIXME: Get rid of this
        # FIXME: port this
        legacy_code = indent(_evaluate_basis_derivatives(ir["evaluate_basis"]), 4)
        return legacy_code


    def evaluate_basis_derivatives_all(self, L, ir, parameters):
        # FIXME: port this
        use_legacy = 1
        if use_legacy:
            return indent(_evaluate_basis_derivatives_all(ir["evaluate_basis"]), 4)

        """
        // Legacy version:
        evaluate_basis_derivatives_all(std::size_t n,
                                       double * values,
                                       const double * x,
                                       const double * coordinate_dofs,
                                       int cell_orientation)
        // Suggestion for new version:
        new_evaluate_basis_derivatives(double * values,
                                       std::size_t order,
                                       std::size_t num_points,
                                       const double * x,
                                       const double * coordinate_dofs,
                                       int cell_orientation,
                                       const ufc::coordinate_mapping * cm)
        """

        # TODO: This is a refactoring step to allow rewriting code
        # generation to use coordinate_mapping in one stage, before
        # making it available as an argument from dolfin in the next stage.
        affine_coordinate_mapping_classname = ir["affine_coordinate_mapping_classname"]

        # Output arguments:
        values = L.Symbol("values")

        # Input arguments:
        #order = L.Symbol("order")
        order = L.Symbol("n")
        x = L.Symbol("x")
        coordinate_dofs = L.Symbol("coordinate_dofs")
        cell_orientation = L.Symbol("cell_orientation")

        # Internal variables:
        #num_points = L.Symbol("num_points")
        num_points = 1  # Always 1 in legacy API
        reference_values = L.Symbol("reference_values")
        X = L.Symbol("X")
        J = L.Symbol("J")
        detJ = L.Symbol("detJ")
        K = L.Symbol("K")
        ip = L.Symbol("ip")

        gdim = ir["geometric_dimension"]
        tdim = ir["topological_dimension"]

        code = [
            # Create local affine coordinate mapping object
            # TODO: Get this as input instead to support non-affine
            L.VariableDecl(affine_coordinate_mapping_classname, "cm"),
            L.ForRange(ip, 0, num_points, index_type=index_type, body=[
                L.ArrayDecl("double", X, (tdim,)),
                L.ArrayDecl("double", J, (gdim*tdim,)),
                L.ArrayDecl("double", detJ, (1,)),
                L.ArrayDecl("double", K, (tdim*gdim,)),
                L.Call("cm.compute_reference_geometry",
                       (X, J, detJ, K, num_points, x, coordinate_dofs, cell_orientation)),
                L.Call("evaluate_reference_basis_derivatives",
                       (reference_values, order, num_points, X)),
                L.Call("transform_reference_basis_derivatives",
                       (values, order, num_points, reference_values, X, J, detJ, K, cell_orientation)),
            ])
        ]
        return code

    def evaluate_dof(self, L, ir, parameters):
        # FIXME: Get rid of this
        # FIXME: port this
        use_legacy = 1
        if use_legacy:
            # Codes generated together
            (legacy_code, evaluate_dofs_code) \
              = evaluate_dof_and_dofs(ir["evaluate_dof"])
        #        print(legacy_code)
        code = _x_evaluate_dof_and_dofs(L, ir["evaluate_dof"])
        # print(L.StatementList(new_code))
        #        return indent(legacy_code, 4)
        return code


    def evaluate_dofs(self, L, ir, parameters):
        """Generate code for evaluate_dofs."""
        """
        - evaluate_dof needs to be split into invert_mapping + evaluate_dof or similar?

          f = M fhat;  nu(f) = nu(M fhat) = nuhat(M^-1 f) = sum_i w_i M^-1 f(x_i)

          // Get fixed set of points on reference element
          num_points = element->num_dof_evaluation_points();
          double X[num_points*tdim];
          element->tabulate_dof_evaluation_points(X);

          // Compute geometry in these points
          domain->compute_geometry(reference_points, num_point, X, J, detJ, K, coordinate_dofs, cell_orientation);

          // Computed by dolfin
          for ip
            fvalues[ip][:] = f.evaluate(point[ip])[:];

          // Finally: nu_j(f) = sum_component sum_ip weights[j][ip][component] fvalues[ip][component]
          element->evaluate_dofs(fdofs, fvalues, J, detJ, K)
        """
        # FIXME: port this, then translate into reference version
        use_legacy = 1
        if use_legacy:
            # Codes generated together
            (evaluate_dof_code, evaluate_dofs_code) \
              = evaluate_dof_and_dofs(ir["evaluate_dof"])
            return indent(evaluate_dofs_code, 4)

    def interpolate_vertex_values(self, L, ir, parameters):
        # legacy_code = indent(interpolate_vertex_values(ir["interpolate_vertex_values"]), 4)
        #        print(legacy_code)

        irdata = ir["interpolate_vertex_values"]
        # Raise error if interpolate_vertex_values is ill-defined
        if not irdata:
            msg = "interpolate_vertex_values is not defined for this element"
            return generate_error(L, msg, parameters["convert_exceptions_to_warnings"])

        # Add code for Jacobian if necessary
        code = []
        gdim = irdata["geometric_dimension"]
        tdim = irdata["topological_dimension"]
        element_cellname = ir["evaluate_basis"]["cellname"]
        if irdata["needs_jacobian"]:
            code += jacobian(L, gdim, tdim, element_cellname)
            code += inverse_jacobian(L, gdim, tdim, element_cellname)
            if irdata["needs_oriented"] and tdim != gdim:
                code += orientation(L)

        # Compute total value dimension for (mixed) element
        total_dim = irdata["physical_value_size"]

        # Generate code for each element
        value_offset = 0
        space_offset = 0
        for data in irdata["element_data"]:
            # Add vertex interpolation for this element
            code += [L.Comment("Evaluate function and change variables")]

            # Extract vertex values for all basis functions
            vertex_values = data["basis_values"]
            value_size = data["physical_value_size"]
            space_dim = data["space_dim"]
            mapping = data["mapping"]

            # Create code for each value dimension:
            for k in range(value_size):
                # Create code for each vertex x_j
                for (j, values_at_vertex) in enumerate(vertex_values):

                    if value_size == 1:
                        values_at_vertex = [values_at_vertex]

                    values = clamp_table_small_numbers(values_at_vertex)

                    # Map basis functions using appropriate mapping
                    # FIXME: sort out all non-affine mappings and make into a function
                    # components = change_of_variables(values_at_vertex, k)

                    if mapping == 'affine':
                        w = values[k]
                    elif mapping == 'contravariant piola':
                        detJ = L.Symbol("detJ")
                        J = L.Symbol("J")
                        w = []
                        for index in range(space_dim):
                            inner = 0.0
                            for p in range(tdim):
                                inner += J[p+k*tdim]*values[p][index]
                            w.append(inner/detJ)
                    elif mapping == 'covariant piola':
                        K = L.Symbol("K")
                        w = []
                        for index in range(space_dim):
                            acc_sum = 0.0
                            for p in range(tdim):
                                acc_sum += K[k+p*gdim]*values[p][index]
                            w.append(acc_sum)
                    elif mapping == 'double covariant piola':
                        K = L.Symbol("K")
                        w = []
                        for index in range(space_dim):
                            acc_sum = 0.0
                            for p in range(tdim):
                                for q in range(tdim):
                                    acc_sum += K[k//tdim + p*gdim]*values[p][q][index]*K[k % tdim + q*gdim]
                            w.append(acc_sum)
                    elif mapping == 'double contravariant piola':
                        J = L.Symbol("J")
                        detJ = L.Symbol("detJ")
                        w = []
                        for index in range(space_dim):
                            acc_sum = 0.0
                            for p in range(tdim):
                                for q in range(tdim):
                                    acc_sum += J[p + (k//tdim)*tdim]*values[p][q][index]*J[q + (k % tdim)*tdim]
                            acc_sum /= (detJ*detJ)
                            w.append(acc_sum)
                    else:
                        raise RuntimeError("Mapping not implemented")

                    # Contract coefficients and basis functions
                    dof_values = L.Symbol("dof_values")
                    dof_list = [dof_values[i + space_offset] for i in range(space_dim)]
                    acc_value = 0.0
                    for p, q in zip(dof_list, w):
                        acc_value += p*q

                    # Assign value to correct vertex
                    index = j * total_dim + (k + value_offset)
                    v_values = L.Symbol("vertex_values")
                    code += [L.Assign(v_values[index], acc_value)]

            # Update offsets for value- and space dimension
            value_offset += data["physical_value_size"]
            space_offset += data["space_dim"]

            #        print(L.StatementList(code))
        return code

    def tabulate_dof_coordinates(self, L, ir, parameters):
        ir = ir["tabulate_dof_coordinates"]

        # Raise error if tabulate_dof_coordinates is ill-defined
        if not ir:
            msg = "tabulate_dof_coordinates is not defined for this element"
            return generate_error(L, msg, parameters["convert_exceptions_to_warnings"])

        # Extract coordinates and cell dimension
        gdim = ir["gdim"]
        tdim = ir["tdim"]
        points = ir["points"]

        # Output argument
        dof_coordinates = L.FlattenedArray(L.Symbol("dof_coordinates"),
                                           dims=(len(points), gdim))

        # Input argument
        coordinate_dofs = L.Symbol("coordinate_dofs")

        # Loop indices
        i = L.Symbol("i")
        k = L.Symbol("k")
        ip = L.Symbol("ip")

        # Basis symbol
        phi = L.Symbol("phi")

        # TODO: Get rid of all places that use affine_weights, assumes affine mesh
        # Create code for evaluating affine coordinate basis functions
        num_scalar_xdofs = tdim + 1
        cg1_basis = affine_weights(tdim)
        phi_values = numpy.asarray([phi_comp for X in points for phi_comp in cg1_basis(X)])
        assert len(phi_values) == len(points) * num_scalar_xdofs

        # TODO: Use precision parameter here
        phi_values = clamp_table_small_numbers(phi_values)

        code = [
            L.Assign(
                dof_coordinates[ip][i],
                sum(phi_values[ip*num_scalar_xdofs + k] * coordinate_dofs[gdim*k + i]
                    for k in range(num_scalar_xdofs))
            )
            for ip in range(len(points))
            for i in range(gdim)
        ]

        # FIXME: This code assumes an affine coordinate field.
        #        To get around that limitation, make this function take another argument
        #            const ufc::coordinate_mapping * cm
        #        and generate code like this:
        """
        index_type X[tdim*num_dofs];
        tabulate_dof_coordinates(X);
        cm->compute_physical_coordinates(x, X, coordinate_dofs);
        """

        return code

    def tabulate_reference_dof_coordinates(self, L, ir, parameters):
        # TODO: Change signature to avoid copy? E.g.
        # virtual const std::vector<double> & tabulate_reference_dof_coordinates() const = 0;
        # See integral::enabled_coefficients for example

        # TODO: ensure points is a numpy array,
        #   get tdim from points.shape[1],
        #   place points in ir directly instead of the subdict
        ir = ir["tabulate_dof_coordinates"]

        # Raise error if tabulate_reference_dof_coordinates is ill-defined
        if not ir:
            msg = "tabulate_reference_dof_coordinates is not defined for this element"
            return generate_error(L, msg, parameters["convert_exceptions_to_warnings"])

        # Extract coordinates and cell dimension
        tdim = ir["tdim"]
        points = ir["points"]

        # Output argument
        reference_dof_coordinates = L.Symbol("reference_dof_coordinates")

        # Reference coordinates
        dof_X = L.Symbol("dof_X")
        dof_X_values = [X[jj] for X in points for jj in range(tdim)]
        decl = L.ArrayDecl("static const double", dof_X,
                           (len(points) * tdim,), values=dof_X_values)
        copy = L.MemCopy(dof_X, reference_dof_coordinates, tdim*len(points))

        code = [decl, copy]
        return code

    def evaluate_reference_basis(self, L, ir, parameters):
        data = ir["evaluate_basis"]
        if isinstance(data, string_types):
            msg = "evaluate_reference_basis: %s" % data
            return generate_error(L, msg, parameters["convert_exceptions_to_warnings"])

        return generate_evaluate_reference_basis(L, data, parameters)

    def evaluate_reference_basis_derivatives(self, L, ir, parameters):
        data = ir["evaluate_basis"]
        if isinstance(data, string_types):
            msg = "evaluate_reference_basis_derivatives: %s" % data
            return generate_error(L, msg, parameters["convert_exceptions_to_warnings"])

        return generate_evaluate_reference_basis_derivatives(L, data, parameters)

    def transform_reference_basis_derivatives(self, L, ir, parameters):
        data = ir["evaluate_basis"]
        if isinstance(data, string_types):
            msg = "transform_reference_basis_derivatives: %s" % data
            return generate_error(L, msg, parameters["convert_exceptions_to_warnings"])

        # Get some known dimensions
        #element_cellname = data["cellname"]
        gdim = data["geometric_dimension"]
        tdim = data["topological_dimension"]
        max_degree = data["max_degree"]
        reference_value_size = data["reference_value_size"]
        physical_value_size = data["physical_value_size"]
        num_dofs = len(data["dofs_data"])

        max_g_d = gdim**max_degree
        max_t_d = tdim**max_degree

        # Output arguments
        values_symbol = L.Symbol("values")

        # Input arguments
        order = L.Symbol("order")
        num_points = L.Symbol("num_points")  # FIXME: Currently assuming 1 point?
        reference_values = L.Symbol("reference_values")
        J = L.Symbol("J")
        detJ = L.Symbol("detJ")
        K = L.Symbol("K")

        # Internal variables
        transform = L.Symbol("transform")

        # Indices, I've tried to use these for a consistent purpose
        ip = L.Symbol("ip") # point
        i = L.Symbol("i")   # physical component
        j = L.Symbol("j")   # reference component
        k = L.Symbol("k")   # order
        r = L.Symbol("r")   # physical derivative number
        s = L.Symbol("s")   # reference derivative number
        d = L.Symbol("d")   # dof

        combinations_code = []
        if max_degree == 0:
            # Don't need combinations
            num_derivatives_t = 1  # TODO: I think this is the right thing to do to make this still work for order=0?
            num_derivatives_g = 1
        elif tdim == gdim:
            num_derivatives_t = L.Symbol("num_derivatives")
            num_derivatives_g = num_derivatives_t
            combinations_code += [
                L.VariableDecl("const " + index_type, num_derivatives_t,
                               L.Call("std::pow", (tdim, order))),
            ]

            # Add array declarations of combinations
            combinations_code_t, combinations_t = _generate_combinations(L, tdim, max_degree, order, num_derivatives_t)
            combinations_code += combinations_code_t
            combinations_g = combinations_t
        else:
            num_derivatives_t = L.Symbol("num_derivatives_t")
            num_derivatives_g = L.Symbol("num_derivatives_g")
            combinations_code += [
                L.VariableDecl("const " + index_type, num_derivatives_t,
                               L.Call("std::pow", (tdim, order))),
                L.VariableDecl("const " + index_type, num_derivatives_g,
                               L.Call("std::pow", (gdim, order))),
            ]
            # Add array declarations of combinations
            combinations_code_t, combinations_t = _generate_combinations(L, tdim, max_degree, order, num_derivatives_t, suffix="_t")
            combinations_code_g, combinations_g = _generate_combinations(L, gdim, max_degree, order, num_derivatives_g, suffix="_g")
            combinations_code += combinations_code_t
            combinations_code += combinations_code_g

        # Define expected dimensions of argument arrays
        J = L.FlattenedArray(J, dims=(num_points, gdim, tdim))
        detJ = L.FlattenedArray(detJ, dims=(num_points,))
        K = L.FlattenedArray(K, dims=(num_points, tdim, gdim))

        values = L.FlattenedArray(values_symbol,
            dims=(num_points, num_dofs, num_derivatives_g, physical_value_size))
        reference_values = L.FlattenedArray(reference_values,
            dims=(num_points, num_dofs, num_derivatives_t, reference_value_size))

        # Generate code to compute the derivative transform matrix
        transform_matrix_code = [
            # Initialize transform matrix to all 1.0
            L.ArrayDecl("double", transform, (max_g_d, max_t_d)),
            L.ForRanges(
                (r, 0, num_derivatives_g),
                (s, 0, num_derivatives_t),
                index_type=index_type,
                body=L.Assign(transform[r, s], 1.0)
            ),
            ]
        if max_degree > 0:
            transform_matrix_code += [
                # Compute transform matrix entries, each a product of K entries
                L.ForRanges(
                    (r, 0, num_derivatives_g),
                    (s, 0, num_derivatives_t),
                    (k, 0, order),
                    index_type=index_type,
                    body=L.AssignMul(transform[r, s],
                                     K[ip, combinations_t[s, k], combinations_g[r, k]])
                ),
            ]

        # Initialize values to 0, will be added to inside loops
        values_init_code = [
            L.MemZero(values_symbol, num_points * num_dofs * num_derivatives_g * physical_value_size),
            ]

        # Make offsets available in generated code
        reference_offsets = L.Symbol("reference_offsets")
        physical_offsets = L.Symbol("physical_offsets")
        dof_attributes_code = [
            L.ArrayDecl("const " + index_type, reference_offsets, (num_dofs,),
                        values=[dof_data["reference_offset"] for dof_data in data["dofs_data"]]),
            L.ArrayDecl("const " + index_type, physical_offsets, (num_dofs,),
                        values=[dof_data["physical_offset"] for dof_data in data["dofs_data"]]),
            ]

        # Build dof lists for each mapping type
        mapping_dofs = defaultdict(list)
        for idof, dof_data in enumerate(data["dofs_data"]):
            mapping_dofs[dof_data["mapping"]].append(idof)

        # Generate code for each mapping type
        d = L.Symbol("d")
        transform_apply_code = []
        for mapping in sorted(mapping_dofs):
            # Get list of dofs using this mapping
            idofs = mapping_dofs[mapping]

            # Select iteration approach over dofs
            if idofs == list(range(idofs[0], idofs[-1]+1)):
                # Contiguous
                dofrange = (d, idofs[0], idofs[-1]+1)
                idof = d
            else:
                # Stored const array of dof indices
                idofs_symbol = L.Symbol("%s_dofs" % mapping.replace(" ", "_"))
                dof_attributes_code += [
                    L.ArrayDecl("const " + index_type, idofs_symbol,
                                (len(idofs),), values=idofs),
                ]
                dofrange = (d, 0, len(idofs))
                idof = idofs_symbol[d]

            # NB! Array access to offsets, these are not Python integers
            reference_offset = reference_offsets[idof]
            physical_offset = physical_offsets[idof]

            # How many components does each basis function with this mapping have?
            # This should be uniform, i.e. there should be only one element in this set:
            num_reference_components, = set(data["dofs_data"][i]["num_components"] for i in idofs)

            M_scale, M_row, num_physical_components = generate_element_mapping(
                mapping, i,
                num_reference_components, tdim, gdim,
                J[ip], detJ[ip], K[ip]
            )

            transform_apply_body = [
                L.AssignAdd(values[ip, idof, r, physical_offset + k],
                            transform[r, s] * reference_values[ip, idof, s, reference_offset + k])
                for k in range(num_physical_components)
            ]

            msg = "Using %s transform to map values back to the physical element." % mapping.replace("piola", "Piola")

            mapped_value = L.Symbol("mapped_value")
            transform_apply_code += [
                L.ForRanges(
                    dofrange,
                    (s, 0, num_derivatives_t),
                    (i, 0, num_physical_components),
                    index_type=index_type, body=[
                        # Unrolled application of mapping to one physical component,
                        # for affine this automatically reduces to
                        #   mapped_value = reference_values[..., reference_offset]
                        L.Comment(msg),
                        L.VariableDecl("const double", mapped_value,
                                       M_scale * sum(M_row[jj] * reference_values[ip, idof, s, reference_offset + jj]
                                                     for jj in range(num_reference_components))),
                        # Apply derivative transformation, for order=0 this reduces to
                        # values[ip,idof,0,physical_offset+i] = transform[0,0]*mapped_value
                        L.Comment("Mapping derivatives back to the physical element"),
                        L.ForRanges(
                            (r, 0, num_derivatives_g),
                            index_type=index_type, body=[
                                L.AssignAdd(values[ip, idof, r, physical_offset + i],
                                            transform[r, s] * mapped_value)
                        ])
                ])
            ]

        # Transform for each point
        point_loop_code = [
            L.ForRange(ip, 0, num_points, index_type=index_type, body=(
                transform_matrix_code
                + transform_apply_code
            ))
        ]

        # Join code
        code = (
            combinations_code
            + values_init_code
            + dof_attributes_code
            + point_loop_code
        )
        return code
