# Copyright (C) 2009-2017 Anders Logg and Martin Sandve Alnæs, Chris Richardson
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

# Note: Much of the code in this file is a direct translation
# from the old implementation in FFC, although some improvements
# have been made to the generated code.

from ffcx.codegeneration.jacobian import inverse_jacobian, jacobian
from ufl.permutation import build_component_numbering

index_type = "int64_t"


def reference_to_physical_map(cellname):
    """Returns a map from reference coordinates to physical element coordinates."""

    if cellname == "interval":
        return lambda x: (1.0 - x[0], x[0])
    elif cellname == "triangle":
        return lambda x: (1.0 - x[0] - x[1], x[0], x[1])
    elif cellname == "tetrahedron":
        return lambda x: (1.0 - x[0] - x[1] - x[2], x[0], x[1], x[2])
    elif cellname == "quadrilateral":
        return lambda x: ((1 - x[0]) * (1 - x[1]), (1 - x[0]) * x[1], x[0] * (1 - x[1]), x[0] * x[1])
    elif cellname == "hexahedron":
        return lambda x: ((1 - x[0]) * (1 - x[1]) * (1 - x[2]), (1 - x[0]) * (1 - x[1]) * x[2], (1 - x[0]) * x[1] * (1 - x[2]), (1 - x[0]) * x[1] * x[2], x[0] * (1 - x[1]) * (1 - x[2]), x[0] * (1 - x[1]) * x[2], x[0] * x[1] * (1 - x[2]), x[0] * x[1] * x[2])  # noqa: E501


def _change_variables(L, mapping, gdim, tdim, offset):
    """Generate code for mapping function values according to 'mapping' and offset.

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

    values = L.Symbol("physical_values")

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
                inner += values[j + offset] * K[i, j]
            w.append(inner * detJ)
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
                inner += values[j + offset] * J[j, i]
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
                        inner += J[j, i] * values[j * tdim + k + offset] * J[k, l]
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
                        inner += K[i, j] * values[j * tdim + k + offset] * K[l, k]
                w.append(inner * detJ * detJ)
        return w

    else:
        raise Exception("The mapping (%s) is not allowed" % mapping)


def _generate_body(L, i, dof, mapping, gdim, tdim, cell_shape, offset=0):
    """Generate code for a single dof."""

    # EnrichedElement is handled by having [None, ..., None] dual basis
    if not dof:
        # Return error code
        return ([L.Return(-1)], 0.0)

    points = list(dof.keys())

    # Generate different code if multiple points. (Otherwise ffcx
    # compile time blows up.)
    if len(points) > 1:
        return ([], -1)
        # FIXME: this currently returns empty code as the line below returns code that fails to compile
        # return _generate_multiple_points_body(L, i, dof, mapping, gdim, tdim, offset)

    # Get weights for mapping reference point to physical
    x = points[0]
    w = reference_to_physical_map(cell_shape)(x)

    # Map point onto physical element: y = F_K(x)
    code = []

    # Map function values to the reference element
    F = _change_variables(L, mapping, gdim, tdim, offset)

    # Simple affine functions deserve special case:
    if len(F) == 1:
        return (code, dof[x][0][0] * F[0])

    # Flatten multi-indices
    (index_map, _) = build_component_numbering([tdim] * len(dof[x][0][1]), ())
    # Take inner product between components and weights
    value = 0.0
    for (w, k) in dof[x]:
        value += w * F[index_map[k]]

    # Return eval code and value
    return (code, value)


def _generate_multiple_points_body(L, i, dof, mapping, gdim, tdim, offset=0):
    """Generate c++ for-loop for multiple points (integral bodies)."""

    result = L.Symbol("result")
    code = [L.Assign(result, 0.0)]
    points = list(dof.keys())
    n = len(points)

    # Get number of tokens per point
    tokens = [dof[x] for x in points]
    len_tokens = [len(t) for t in tokens]
    assert all(x == len_tokens[0] for x in len_tokens)
    len_tokens = len_tokens[0]

    # Declare points
    #    points = format["list"]([format["list"](x) for x in points])
    X_i = L.Symbol("X_%d" % i)
    code += [L.ArrayDecl("double", X_i, [n, tdim], points)]

    # Declare components
    components = [[c[0] for (w, c) in token] for token in tokens]
    #    components = format["list"]([format["list"](c) for c in components])
    D_i = L.Symbol("D_%d" % i)
    code += [L.ArrayDecl("int", D_i, [n, len_tokens], components)]

    # Declare weights
    weights = [[w for (w, c) in token] for token in tokens]
    #    weights = format["list"]([format["list"](w) for w in weights])
    W_i = L.Symbol("W_%d" % i)
    code += [L.ArrayDecl("double", W_i, [n, len_tokens], weights)]

    # Declare copy variable:
    copy_i = L.Symbol("copy_%d" % i)
    code += [L.ArrayDecl("double", copy_i, tdim)]

    # Add loop over points
    code += [L.Comment("Loop over points")]

    # Map the points from the reference onto the physical element
    r = L.Symbol("r")
    w0 = L.Symbol("w0")
    w1 = L.Symbol("w1")
    w2 = L.Symbol("w2")
    w3 = L.Symbol("w3")
    y = L.Symbol("y")

    coordinate_dofs = L.Symbol("coordinate_dofs")

    if tdim == 1:
        lines_r = [
            L.Comment("Evaluate basis functions for affine mapping"),
            L.VariableDecl("const double", w0, 1 - X_i[r][0]),
            L.VariableDecl("const double", w1, X_i[r][0]),
            L.Comment("Compute affine mapping y = F(X)"),
            L.ArrayDecl("double", y, [gdim])
        ]
        for j in range(gdim):
            lines_r += [L.Assign(y[j], w0 * coordinate_dofs[j] + w1 * coordinate_dofs[j + gdim])]
    elif tdim == 2:
        lines_r = [
            L.Comment("Evaluate basis functions for affine mapping"),
            L.VariableDecl("const double", w0, 1 - X_i[r][0] - X_i[r][1]),
            L.VariableDecl("const double", w1, X_i[r][0]),
            L.VariableDecl("const double", w2, X_i[r][1]),
            L.Comment("Compute affine mapping y = F(X)"),
            L.ArrayDecl("double", y, [gdim])
        ]
        for j in range(gdim):
            lines_r += [
                L.Assign(y[j], w0 * coordinate_dofs[j]
                         + w1 * coordinate_dofs[j + gdim]
                         + w2 * coordinate_dofs[j + 2 * gdim])
            ]
    elif tdim == 3:
        lines_r = [
            L.Comment("Evaluate basis functions for affine mapping"),
            L.VariableDecl("const double", w0, 1 - X_i[r][0] - X_i[r][1] - X_i[r][2]),
            L.VariableDecl("const double", w1, X_i[r][0]),
            L.VariableDecl("const double", w2, X_i[r][1]),
            L.VariableDecl("const double", w3, X_i[r][2]),
            L.Comment("Compute affine mapping y = F(X)"),
            L.ArrayDecl("double", y, [gdim])
        ]
        for j in range(gdim):
            lines_r += [
                L.Assign(y[j], w0 * coordinate_dofs[j]
                         + w1 * coordinate_dofs[j + gdim]
                         + w2 * coordinate_dofs[j + 2 * gdim]
                         + w3 * coordinate_dofs[j + 3 * gdim])
            ]

    # Evaluate function at physical point
    lines_r += [L.Comment("Evaluate function at physical point")]
    # vals = L.Symbol("physical_values")
    # c = L.Symbol("c")
    lines_r += [L.Comment("Replace this: f.evaluate(vals, y, c)")]

    # Map function values to the reference element
    lines_r += [L.Comment("Map function to reference element")]
    F = _change_variables(L, mapping, gdim, tdim, offset)
    lines_r += [L.Assign(copy_i[k], F_k) for (k, F_k) in enumerate(F)]

    # Add loop over directional components
    lines_r += [L.Comment("Loop over directions")]

    s = L.Symbol("s")
    lines_r += [
        L.ForRange(
            s,
            0,
            len_tokens,
            index_type=index_type,
            body=[L.AssignAdd(result, copy_i[D_i[r, s]] * W_i[r, s])])
    ]

    # Generate loop over r and add to code.
    code += [L.ForRange(r, 0, n, index_type=index_type, body=lines_r)]
    return (code, result)


def generate_transform_values(L, ir):
    """Generate code for transform_values.
    Transforms values in physical space into reference space.
    These values represent evaluation of the function at dof
    points (only valid for point evaluation dofs).

    """

    gdim = ir.geometric_dimension
    tdim = ir.topological_dimension
    cell_shape = ir.cell_shape

    # Enriched element, no dofs defined
    if not any(ir.dofs):
        code = []
    else:
        code = []
        # Check whether Jacobians are necessary.
        needs_inverse_jacobian = any(["contravariant piola" in m for m in ir.mappings])
        needs_jacobian = any(["covariant piola" in m for m in ir.mappings])

        # Intermediate variable needed for multiple point dofs
        needs_temporary = any(dof is not None and len(dof) > 1 for dof in ir.dofs)
        if needs_temporary:
            result = L.Symbol("result")
            code += [L.VariableDecl("double", result)]

        if needs_jacobian or needs_inverse_jacobian:
            code += jacobian(L, gdim, tdim, cell_shape)

        if needs_inverse_jacobian:
            code += inverse_jacobian(L, gdim, tdim, cell_shape)

    # Extract variables
    mappings = ir.mappings
    offsets = ir.physical_offsets

    # Generate bodies for each degree of freedom
    values = L.Symbol("reference_values")
    value_size = ir.physical_value_size
    for (i, dof) in enumerate(ir.dofs):
        c, r = _generate_body(L, i, dof, mappings[i], gdim, tdim, cell_shape,
                              offsets[i] + i * value_size)
        code += c
        code += [L.Assign(values[i], r)]

    code += [L.Return(0)]
    return code
