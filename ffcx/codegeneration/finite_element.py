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
from ffcx.codegeneration.utils import (generate_return_int_switch,
                                       generate_return_new_switch,
                                       apply_transformations_to_data)

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

    d["value_dimension"] = value_dimension(L, ir.value_shape)
    d["reference_value_dimension"] = reference_value_dimension(L, ir.reference_value_shape)

    d["apply_dof_transformation"] = L.StatementList(
        apply_dof_transformation(L, ir, parameters))
    d["apply_dof_transformation_to_scalar"] = L.StatementList(
        apply_dof_transformation(L, ir, parameters, dtype="ufc_scalar_t"))
    d["apply_inverse_transpose_dof_transformation"] = L.StatementList(
        apply_dof_transformation(L, ir, parameters, inverse=True, transpose=True))
    d["apply_inverse_transpose_dof_transformation_to_scalar"] = L.StatementList(
        apply_dof_transformation(L, ir, parameters, dtype="ufc_scalar_t", inverse=True, transpose=True))

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

    # Define space dimension for the element
    if parameters.get("sycl_defines", False):
        elm_defininiton = "space_dimension_" + ir.name
        define_dimension = "\n#define " + elm_defininiton + " " + str(ir.space_dimension) + "\n"
        implementation = define_dimension + implementation

    return declaration, implementation
