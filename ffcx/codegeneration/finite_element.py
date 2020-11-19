# Copyright (C) 2009-2017 Anders Logg and Martin Sandve Alnæs
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

# Note: Much of the code in this file is a direct translation
# from the old implementation in FFC, although some improvements
# have been made to the generated code.

from collections import defaultdict
import logging

import numpy
import ffcx.codegeneration.finite_element_template as ufc_finite_element
import ufl
from ffcx.codegeneration.utils import (generate_return_int_switch,
                                       generate_return_new_switch)

logger = logging.getLogger("ffcx")

index_type = "int"


def _generate_combinations(L, tdim, max_degree, order, num_derivatives, suffix=""):
    max_num_derivatives = tdim**max_degree
    combinations = L.Symbol("combinations" + suffix)

    # This precomputes the combinations for each order and stores in code as table
    # Python equivalent precomputed for each valid order:
    combinations_shape = (max_degree, max_num_derivatives, max_degree)
    all_combinations = numpy.zeros(combinations_shape, dtype=int)
    for q in range(1, max_degree + 1):
        for row in range(1, max_num_derivatives):
            for num in range(0, row):
                for col in range(q - 1, -1, -1):
                    if all_combinations[q - 1][row][col] > tdim - 2:
                        all_combinations[q - 1][row][col] = 0
                    else:
                        all_combinations[q - 1][row][col] += 1
                        break
    code = [
        L.Comment("Precomputed combinations"),
        L.ArrayDecl(
            "const " + index_type, combinations, combinations_shape, values=all_combinations),
    ]
    # Select the right order for further access
    combinations = combinations[order - 1]

    return code, combinations


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
        M_scale = 1.0
        M_row = [K[jj, i] for jj in range(tdim)]
    elif mapping == "double covariant piola":
        assert num_reference_components == tdim**2
        num_physical_components = gdim**2
        # g_il = K_ji G_jk K_kl = K_ji K_kl G_jk
        i0 = i // tdim  # i in the line above
        i1 = i % tdim  # l ...
        M_scale = 1.0
        M_row = [K[jj, i0] * K[kk, i1] for jj in range(tdim) for kk in range(tdim)]
    elif mapping == "double contravariant piola":
        assert num_reference_components == tdim**2
        num_physical_components = gdim**2
        # g_il = (det J)^(-2) Jij G_jk Jlk = (det J)^(-2) Jij Jlk G_jk
        i0 = i // tdim  # i in the line above
        i1 = i % tdim  # l ...
        M_scale = 1.0 / (detJ * detJ)
        M_row = [J[i0, jj] * J[i1, kk] for jj in range(tdim) for kk in range(tdim)]
    else:
        raise RuntimeError("Unknown mapping: %s" % mapping)

    return M_scale, M_row, num_physical_components


def value_dimension(L, value_shape):
    return generate_return_int_switch(L, "i", value_shape, 1)


def reference_value_dimension(L, reference_value_shape):
    return generate_return_int_switch(L, "i", reference_value_shape, 1)


def sub_element_declaration(L, ir):
    classnames = set(ir.create_sub_element)
    code = ""
    for name in classnames:
        code += f"ufc_finite_element* create_{name}(void);\n"
    return code


def create_sub_element(L, ir):
    classnames = ir.create_sub_element
    return generate_return_new_switch(L, "i", classnames)


def transform_values(L, ir, parameters):
    """Generate code for transform_values."""
    return [L.Return(-1)]  # generate_transform_values(L, ir.evaluate_dof)


def tabulate_reference_dof_coordinates(L, ir, parameters):
    # TODO: ensure points is a numpy array,
    #   get tdim from points.shape[1],
    #   place points in ir directly instead of the subdict
    ir = ir.tabulate_dof_coordinates

    # Raise error if tabulate_reference_dof_coordinates is ill-defined
    if not ir:
        return [L.Return(-1)]

    # Extract coordinates and cell dimension
    tdim = ir.tdim
    points = ir.points

    # Output argument
    reference_dof_coordinates = L.Symbol("reference_dof_coordinates")

    # Reference coordinates
    dof_X = L.Symbol("dof_X")
    dof_X_values = [X[jj] for X in points for jj in range(tdim)]
    decl = L.ArrayDecl("static const double", dof_X, (len(points) * tdim, ), values=dof_X_values)
    copy = L.MemCopy(dof_X, reference_dof_coordinates, tdim * len(points), "double")
    ret = L.Return(0)
    return [decl, copy, ret]


def entity_reflection(L, i, cell_shape):
    """Returns the bool that says whether or not an entity has been reflected."""
    cell_info = L.Symbol("cell_permutation")
    if cell_shape in ["triangle", "quadrilateral"]:
        num_faces = 0
        face_bitsize = 1
        assert i[0] == 1
    if cell_shape == "tetrahedron":
        num_faces = 4
        face_bitsize = 3
    if cell_shape == "hexahedron":
        num_faces = 6
        face_bitsize = 3
    if i[0] == 1:
        return L.NE(L.BitwiseAnd(cell_info, L.BitShiftL(1, face_bitsize * num_faces + i[1])), 0)
    elif i[0] == 2:
        return L.NE(L.BitwiseAnd(cell_info, L.BitShiftL(1, face_bitsize * i[1])), 0)
    return L.LiteralBool(False)


def entity_rotations(L, i, cell_shape):
    """Returns number of times an entity has been rotates."""
    cell_info = L.Symbol("cell_permutation")
    assert cell_shape in ["tetrahedron", "hexahedron"]
    assert i[0] == 2
    return L.BitwiseAnd(L.BitShiftR(cell_info, 3 * i[1]), 7)


def transform_reference_basis_derivatives(L, ir, parameters):
    # Get some known dimensions
    # element_cellname = data["cellname"]
    gdim = ir.geometric_dimension
    tdim = ir.topological_dimension
    max_degree = ir.degree
    reference_value_size = ufl.product(ir.reference_value_shape)
    physical_value_size = ufl.product(ir.value_shape)
    num_dofs = ir.space_dimension

    max_g_d = gdim**max_degree
    max_t_d = tdim**max_degree

    # Output arguments
    values_symbol = L.Symbol("values")

    # Input arguments
    order = L.Symbol("order")
    # FIXME: Currently assuming 1 point?
    num_points = L.Symbol("num_points")
    reference_values = L.Symbol("reference_values")
    J = L.Symbol("J")
    detJ = L.Symbol("detJ")
    K = L.Symbol("K")

    # Internal variables
    transform = L.Symbol("transform")

    # Indices, I've tried to use these for a consistent purpose
    ip = L.Symbol("ip")  # point
    i = L.Symbol("i")  # physical component
    k = L.Symbol("k")  # order
    r = L.Symbol("r")  # physical derivative number
    s = L.Symbol("s")  # reference derivative number
    d = L.Symbol("d")  # dof

    iz = L.Symbol("l")  # zeroing arrays

    combinations_code = []
    if max_degree == 0:
        # Don't need combinations
        # TODO: I think this is the right thing to do to make this still work for order=0?
        num_derivatives_t = 1
        num_derivatives_g = 1
    elif tdim == gdim:
        num_derivatives_t = L.Symbol("num_derivatives")
        num_derivatives_g = num_derivatives_t
        combinations_code += [
            L.VariableDecl("const " + index_type, num_derivatives_t, L.Call("pow", (tdim, order))),
        ]

        # Add array declarations of combinations
        combinations_code_t, combinations_t = _generate_combinations(L, tdim, max_degree, order,
                                                                     num_derivatives_t)
        combinations_code += combinations_code_t
        combinations_g = combinations_t
    else:
        num_derivatives_t = L.Symbol("num_derivatives_t")
        num_derivatives_g = L.Symbol("num_derivatives_g")
        combinations_code += [
            L.VariableDecl("const " + index_type, num_derivatives_t, L.Call("pow", (tdim, order))),
            L.VariableDecl("const " + index_type, num_derivatives_g, L.Call("pow", (gdim, order))),
        ]
        # Add array declarations of combinations
        combinations_code_t, combinations_t = _generate_combinations(
            L, tdim, max_degree, order, num_derivatives_t, suffix="_t")
        combinations_code_g, combinations_g = _generate_combinations(
            L, gdim, max_degree, order, num_derivatives_g, suffix="_g")
        combinations_code += combinations_code_t
        combinations_code += combinations_code_g

    # Define expected dimensions of argument arrays
    J = L.FlattenedArray(J, dims=(num_points, gdim, tdim))
    detJ = L.FlattenedArray(detJ, dims=(num_points, ))
    K = L.FlattenedArray(K, dims=(num_points, tdim, gdim))

    values = L.FlattenedArray(
        values_symbol, dims=(num_points, num_dofs, num_derivatives_g, physical_value_size))
    reference_values = L.FlattenedArray(
        reference_values, dims=(num_points, num_dofs, num_derivatives_t, reference_value_size))

    # Generate code to compute the derivative transform matrix
    transform_matrix_code = [
        # Initialize transform matrix to all 1.0
        L.ArrayDecl("double", transform, (max_g_d, max_t_d)),
        L.ForRanges((r, 0, num_derivatives_g), (s, 0, num_derivatives_t),
                    index_type=index_type,
                    body=L.Assign(transform[r, s], 1.0)),
    ]
    if max_degree > 0:
        transform_matrix_code += [
            # Compute transform matrix entries, each a product of K entries
            L.ForRanges((r, 0, num_derivatives_g), (s, 0, num_derivatives_t), (k, 0, order),
                        index_type=index_type,
                        body=L.AssignMul(transform[r, s],
                                         K[ip, combinations_t[s, k], combinations_g[r, k]])),
        ]

    # Initialize values to 0, will be added to inside loops
    values_init_code = [
        L.ForRange(
            iz,
            0,
            num_points * num_dofs * num_derivatives_g * physical_value_size,
            index_type=index_type,
            body=L.Assign(values_symbol[iz], 0.0)),
    ]

    # Make offsets available in generated code
    reference_offsets = L.Symbol("reference_offsets")
    physical_offsets = L.Symbol("physical_offsets")
    dof_attributes_code = [
        L.ArrayDecl(
            "const " + index_type,
            reference_offsets, (num_dofs, ),
            values=ir.reference_offsets),
        L.ArrayDecl(
            "const " + index_type,
            physical_offsets, (num_dofs, ),
            values=ir.physical_offsets)
    ]

    # Build dof lists for each mapping type
    mapping_dofs = defaultdict(list)
    for idof, mapping in enumerate(ir.dof_mappings):
        mapping_dofs[mapping].append(idof)

    # Generate code for each mapping type
    d = L.Symbol("d")
    transform_apply_code = []
    for mapping in sorted(mapping_dofs):
        # Get list of dofs using this mapping
        idofs = mapping_dofs[mapping]

        # Select iteration approach over dofs
        if idofs == list(range(idofs[0], idofs[-1] + 1)):
            # Contiguous
            dofrange = (d, idofs[0], idofs[-1] + 1)
            idof = d
        else:
            # Stored const array of dof indices
            idofs_symbol = L.Symbol("%s_dofs" % mapping.replace(" ", "_"))
            dof_attributes_code += [
                L.ArrayDecl("const " + index_type, idofs_symbol, (len(idofs), ), values=idofs),
            ]
            dofrange = (d, 0, len(idofs))
            idof = idofs_symbol[d]

        # NB! Array access to offsets, these are not Python integers
        reference_offset = reference_offsets[idof]
        physical_offset = physical_offsets[idof]

        # How many components does each basis function with this mapping have?
        # This should be uniform, i.e. there should be only one element in this set:
        num_reference_components = ir.num_reference_components

        M_scale, M_row, num_physical_components = generate_element_mapping(
            mapping, i, num_reference_components, tdim, gdim, J[ip], detJ[ip], K[ip])

        #            transform_apply_body = [
        #                L.AssignAdd(values[ip, idof, r, physical_offset + k],
        #                            transform[r, s] * reference_values[ip, idof, s, reference_offset + k])
        #                for k in range(num_physical_components)
        #            ]

        msg = "Using %s transform to map values back to the physical element." % mapping.replace(
            "piola", "Piola")

        mapped_value = L.Symbol("mapped_value")

        transform_apply_code += [
            L.ForRanges(
                dofrange,
                (s, 0, num_derivatives_t),
                (i, 0, num_physical_components),
                index_type=index_type,
                body=[
                    # Unrolled application of mapping to one physical component,
                    # for affine this automatically reduces to
                    #   mapped_value = reference_values[..., reference_offset]
                    L.Comment(msg),
                    L.VariableDecl(
                        "const double", mapped_value,
                        M_scale * sum(
                            M_row[jj] * reference_values[ip, idof, s, reference_offset + jj]
                            for jj in range(num_reference_components))),
                    # Apply derivative transformation, for order=0 this reduces to
                    # values[ip,idof,0,physical_offset+i] = transform[0,0]*mapped_value
                    L.Comment("Mapping derivatives back to the physical element"),
                    L.ForRanges((r, 0, num_derivatives_g),
                                index_type=index_type,
                                body=[
                                    L.AssignAdd(values[ip, idof, r, physical_offset + i],
                                                transform[r, s] * mapped_value)])
                ])
        ]

    base_perms = ir.base_permutations

    if ir.cell_shape == "interval":
        entities = {}
    elif ir.cell_shape == "triangle":
        entities = {1: 3}
    elif ir.cell_shape == "quadrilateral":
        entities = {1: 4}
    elif ir.cell_shape == "tetrahedron":
        entities = {1: 6, 2: 4}
        face_rotation_order = 3
    elif ir.cell_shape == "hexahedron":
        entities = {1: 12, 2: 6}
        face_rotation_order = 4
    else:
        raise NotImplementedError

    perm_n = 0
    perm_data = []
    if 1 in entities:
        for edge in range(entities[1]):
            perm_data.append((
                entity_reflection(L, (1, edge), ir.cell_shape),
                None,
                base_perms[perm_n]
            ))
            perm_n += 1
    if 2 in entities:
        for face in range(entities[2]):
            for rot in range(1, face_rotation_order):
                perm_data.append((
                    entity_rotations(L, (2, face), ir.cell_shape),
                    rot,
                    numpy.linalg.matrix_power(base_perms[perm_n], rot)
                ))
            perm_n += 1
            perm_data.append((
                entity_reflection(L, (2, face), ir.cell_shape),
                None,
                base_perms[perm_n]
            ))
            perm_n += 1

    assert perm_n == len(base_perms)

    # Apply entity permutations
    apply_permutations = []
    temporary_variables = 0
    for entity_perm, value, perm in perm_data:
        body = []

        # Use temporary variables t0, t1, ... to store current data
        temps = {}
        for index, row in enumerate(perm):
            if not numpy.allclose(row, [1 if i == index else 0 for i, j in enumerate(row)]):
                for dof, w in enumerate(row):
                    if not numpy.isclose(w, 0) and dof not in temps:
                        temps[dof] = L.Symbol("t" + str(len(temps)))
                body.append(L.Assign(values[ip, dof, r, physical_offsets[dof] + i],
                                     sum(w * temps[dof] for dof, w in enumerate(row) if not numpy.isclose(w, 0))))
        temporary_variables = max(temporary_variables, len(temps))

        # If changes would be made, continue to next entity
        if len(body) == 0:
            continue

        if value is None:
            condition = entity_perm
        else:
            condition = L.EQ(entity_perm, value)

        body = [L.Assign(t, values[ip, dof, r, physical_offsets[dof] + i]) for dof, t in temps.items()] + body
        apply_permutations.append(L.If(condition,
                                       L.ForRanges((s, 0, num_derivatives_t), (i, 0, num_physical_components),
                                                   (r, 0, num_derivatives_g), index_type=index_type, body=body)))

    if len(apply_permutations) > 0:
        apply_permutations = [L.VariableDecl("double", L.Symbol("t" + str(i)), 0)
                              for i in range(temporary_variables)] + apply_permutations

    # Transform for each point
    point_loop_code = [
        L.ForRange(
            ip,
            0,
            num_points,
            index_type=index_type,
            body=(transform_matrix_code + transform_apply_code + apply_permutations))
    ]

    # Join code
    code = (combinations_code + values_init_code + dof_attributes_code + point_loop_code
            + [L.Comment(msg), L.Return(0)])

    return code


def generator(ir, parameters):
    """Generate UFC code for a finite element."""

    logger.info("Generating code for finite element:")
    logger.info(f"--- family: {ir.family}")
    logger.info(f"--- degree: {ir.degree}")
    logger.info(f"--- value shape: {ir.value_shape}")
    logger.info(f"--- name: {ir.name}")

    d = {}
    d["factory_name"] = ir.name
    d["signature"] = f"\"{ir.signature}\""
    d["geometric_dimension"] = ir.geometric_dimension
    d["topological_dimension"] = ir.topological_dimension
    d["cell_shape"] = ir.cell_shape
    d["space_dimension"] = ir.space_dimension
    d["value_rank"] = len(ir.value_shape)
    d["value_size"] = ufl.product(ir.value_shape)
    d["reference_value_rank"] = len(ir.reference_value_shape)
    d["reference_value_size"] = ufl.product(ir.reference_value_shape)
    d["degree"] = ir.degree
    d["family"] = f"\"{ir.family}\""
    d["num_sub_elements"] = ir.num_sub_elements
    d["block_size"] = ir.block_size

    import ffcx.codegeneration.C.cnodes as L

    d["value_dimension"] = value_dimension(L, ir.value_shape)
    d["reference_value_dimension"] = reference_value_dimension(L, ir.reference_value_shape)

    statements = transform_reference_basis_derivatives(L, ir, parameters)
    d["transform_reference_basis_derivatives"] = L.StatementList(statements)

    statements = transform_values(L, ir, parameters)
    d["transform_values"] = L.StatementList(statements)

    statements = tabulate_reference_dof_coordinates(L, ir, parameters)
    d["tabulate_reference_dof_coordinates"] = L.StatementList(statements)

    statements = create_sub_element(L, ir)
    d["sub_element_declaration"] = sub_element_declaration(L, ir)
    d["create_sub_element"] = statements

    # Check that no keys are redundant or have been missed
    from string import Formatter
    fieldnames = [
        fname for _, fname, _, _ in Formatter().parse(ufc_finite_element.factory) if fname
    ]
    assert set(fieldnames) == set(
        d.keys()), "Mismatch between keys in template and in formattting dict"

    # Format implementation code
    implementation = ufc_finite_element.factory.format_map(d)

    # Format declaration
    declaration = ufc_finite_element.declaration.format(factory_name=ir.name)

    return declaration, implementation
