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
                                       generate_return_new_switch,
                                       apply_permutations_to_data)
from ffcx.basix_interface import MappingType, mapping_to_str

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
    if mapping == MappingType.identity:
        assert num_reference_components == 1
        num_physical_components = 1
        M_scale = 1
        M_row = [1]  # M_row[0] == 1
    elif mapping == MappingType.contravariantPiola:
        assert num_reference_components == tdim
        num_physical_components = gdim
        M_scale = 1.0 / detJ
        M_row = [J[i, jj] for jj in range(tdim)]
    elif mapping == MappingType.covariantPiola:
        assert num_reference_components == tdim
        num_physical_components = gdim
        M_scale = 1.0
        M_row = [K[jj, i] for jj in range(tdim)]
    elif mapping == MappingType.doubleCovariantPiola:
        assert num_reference_components == tdim**2
        num_physical_components = gdim**2
        # g_il = K_ji G_jk K_kl = K_ji K_kl G_jk
        i0 = i // tdim  # i in the line above
        i1 = i % tdim  # l ...
        M_scale = 1.0
        M_row = [K[jj, i0] * K[kk, i1] for jj in range(tdim) for kk in range(tdim)]
    elif mapping == MappingType.doubleContravariantPiola:
        assert num_reference_components == tdim**2
        num_physical_components = gdim**2
        # g_il = (det J)^(-2) Jij G_jk Jlk = (det J)^(-2) Jij Jlk G_jk
        i0 = i // tdim  # i in the line above
        i1 = i % tdim  # l ...
        M_scale = 1.0 / (detJ * detJ)
        M_row = [J[i0, jj] * J[i1, kk] for jj in range(tdim) for kk in range(tdim)]
    else:
        raise RuntimeError("Unknown mapping")

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


def apply_dof_transformation(L, ir, parameters, reverse=False, dtype="double"):
    """Write function that applies the DOF tranformations/permutations to some data."""
    data = L.Symbol("data")
    block = L.Symbol("block")
    block_size = L.Symbol("dim")

    apply_permutations = apply_permutations_to_data(
        L, ir.base_permutations, ir.cell_shape, data, reverse=reverse,
        indices=lambda dof: dof * block_size + block, ranges=[(block, 0, block_size)],
        dtype=dtype)
    return apply_permutations + [L.Return(0)]


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
    d["needs_permutation_data"] = ir.needs_permutation_data
    d["interpolation_is_identity"] = ir.interpolation_is_identity

    import ffcx.codegeneration.C.cnodes as L

    d["value_dimension"] = value_dimension(L, ir.value_shape)
    d["reference_value_dimension"] = reference_value_dimension(L, ir.reference_value_shape)

    statements = apply_dof_transformation(L, ir, parameters)
    d["apply_dof_transformation"] = L.StatementList(statements)

    statements = apply_dof_transformation(L, ir, parameters, dtype="ufc_scalar_t")
    d["apply_dof_transformation_to_scalar"] = L.StatementList(statements)

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
