# -*- coding: utf-8 -*-
# Copyright (C) 2009-2017 Anders Logg and Martin Sandve Alnæs
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

# Note: Much of the code in this file is a direct translation
# from the old implementation in FFC, although some improvements
# have been made to the generated code.

from collections import defaultdict

import ffc.codegeneration.finite_element_template as ufc_finite_element
import ufl
from ffc.codegeneration.evalderivs import (_generate_combinations,
                                           generate_evaluate_reference_basis_derivatives)
from ffc.codegeneration.evaluatebasis import generate_evaluate_reference_basis
from ffc.codegeneration.evaluatedof import generate_transform_values
from ffc.codegeneration.utils import (generate_return_int_switch,
                                      generate_return_new_switch)

index_type = "int64_t"


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
        code += "ufc_finite_element* create_{name}(void);\n".format(name=name)
    return code


def create_sub_element(L, ir):
    classnames = ir.create_sub_element
    return generate_return_new_switch(L, "i", classnames)


def transform_values(L, ir, parameters):
    """Generate code for transform_values()"""
    return generate_transform_values(L, ir.evaluate_dof)


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


def evaluate_reference_basis(L, ir, parameters):
    data = ir.evaluate_basis
    if isinstance(data, str):
        # Function has not been requested
        msg = "evaluate_reference_basis: {}".format(data)
        return [L.Comment(msg), L.Return(-1)]

    return generate_evaluate_reference_basis(L, data, parameters)


def evaluate_reference_basis_derivatives(L, ir, parameters):
    data = ir.evaluate_basis
    if isinstance(data, str):
        # Function has not been requested
        msg = "evaluate_reference_basis_derivatives: {}".format(data)
        return [L.Comment(msg), L.Return(-1)]

    return generate_evaluate_reference_basis_derivatives(L, data, ir.classname, parameters)


def transform_reference_basis_derivatives(L, ir, parameters):
    data = ir.evaluate_basis
    if isinstance(data, str):
        # Function has not been requested
        msg = "transform_reference_basis_derivatives: {}".format(data)
        return [L.Comment(msg), L.Return(-1)]

    # Get some known dimensions
    # element_cellname = data["cellname"]
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
            values=[dof_data["reference_offset"] for dof_data in data["dofs_data"]]),
        L.ArrayDecl(
            "const " + index_type,
            physical_offsets, (num_dofs, ),
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
        num_reference_components, = set(data["dofs_data"][i]["num_components"] for i in idofs)

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

    # Transform for each point
    point_loop_code = [
        L.ForRange(
            ip,
            0,
            num_points,
            index_type=index_type,
            body=(transform_matrix_code + transform_apply_code))
    ]

    # Join code
    code = (combinations_code + values_init_code + dof_attributes_code + point_loop_code
            + [L.Comment(msg), L.Return(0)])
    return code


def generator(ir, parameters):
    """Generate UFC code for a finite element"""
    d = {}
    d["factory_name"] = ir.classname
    d["signature"] = "\"{}\"".format(ir.signature)
    d["geometric_dimension"] = ir.geometric_dimension
    d["topological_dimension"] = ir.topological_dimension
    d["cell_shape"] = ir.cell_shape
    d["space_dimension"] = ir.space_dimension
    d["value_rank"] = len(ir.value_shape)
    d["value_size"] = ufl.product(ir.value_shape)
    d["reference_value_rank"] = len(ir.reference_value_shape)
    d["reference_value_size"] = ufl.product(ir.reference_value_shape)
    d["degree"] = ir.degree
    d["family"] = "\"{}\"".format(ir.family)
    d["num_sub_elements"] = ir.num_sub_elements

    import ffc.codegeneration.C.cnodes as L

    d["value_dimension"] = value_dimension(L, ir.value_shape)
    d["reference_value_dimension"] = reference_value_dimension(L, ir.reference_value_shape)

    statements = evaluate_reference_basis(L, ir, parameters)
    d["evaluate_reference_basis"] = L.StatementList(statements)

    statements = evaluate_reference_basis_derivatives(L, ir, parameters)
    d["evaluate_reference_basis_derivatives"] = L.StatementList(statements)

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
    declaration = ufc_finite_element.declaration.format(factory_name=ir.classname)

    return declaration, implementation
