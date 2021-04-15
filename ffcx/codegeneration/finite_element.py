# Copyright (C) 2009-2017 Anders Logg and Martin Sandve Alnæs
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

# Note: Much of the code in this file is a direct translation
# from the old implementation in FFC, although some improvements
# have been made to the generated code.

import logging

import numpy
import ffcx.codegeneration.finite_element_template as ufc_finite_element
import ufl
from ffcx.codegeneration.utils import (apply_transformations_to_data)

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


def apply_dof_transformation(L, ir, parameters, inverse=False, transpose=False, dtype="double"):
    """Write function that applies the DOF transformations to some data."""
    data = L.Symbol("data")
    block = L.Symbol("block")
    block_size = L.Symbol("dim")

    apply_transformations = apply_transformations_to_data(
        L, ir.base_transformations, ir.cell_shape, data, inverse=inverse, transpose=transpose,
        indices=lambda dof: dof * block_size + block, ranges=[(block, 0, block_size)],
        dtype=dtype)
    return apply_transformations + [L.Return(0)]


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
    d["needs_transformation_data"] = ir.needs_transformation_data
    d["interpolation_is_identity"] = ir.interpolation_is_identity

    import ffcx.codegeneration.C.cnodes as L

    if len(ir.value_shape) > 0:
        d["value_shape"] = f"value_shape_{ir.name}"
        d["value_shape_init"] = L.ArrayDecl(
            "int", f"value_shape_{ir.name}", values=ir.value_shape, sizes=len(ir.value_shape))
    else:
        d["value_shape"] = "NULL"
        d["value_shape_init"] = ""

    if len(ir.value_shape) > 0:
        d["reference_value_shape"] = f"reference_value_shape_{ir.name}"
        d["reference_value_shape_init"] = L.ArrayDecl(
            "int", f"reference_value_shape_{ir.name}",
            values=ir.reference_value_shape, sizes=len(ir.reference_value_shape))
    else:
        d["reference_value_shape"] = "NULL"
        d["reference_value_shape_init"] = ""

    d["apply_dof_transformation"] = L.StatementList(
        apply_dof_transformation(L, ir, parameters))
    d["apply_dof_transformation_to_scalar"] = L.StatementList(
        apply_dof_transformation(L, ir, parameters, dtype="ufc_scalar_t"))
    d["apply_inverse_transpose_dof_transformation"] = L.StatementList(
        apply_dof_transformation(L, ir, parameters, inverse=True, transpose=True))
    d["apply_inverse_transpose_dof_transformation_to_scalar"] = L.StatementList(
        apply_dof_transformation(L, ir, parameters, dtype="ufc_scalar_t", inverse=True, transpose=True))

    if len(ir.sub_elements) > 0:
        d["sub_elements"] = f"sub_elements_{ir.name}"
        d["sub_elements_init"] = L.ArrayDecl(
            "ufc_finite_element*", f"sub_elements_{ir.name}",
            values=[L.AddressOf(L.Symbol(el)) for el in ir.sub_elements], sizes=len(ir.sub_elements))
    else:
        d["sub_elements"] = "NULL"
        d["sub_elements_init"] = ""

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
