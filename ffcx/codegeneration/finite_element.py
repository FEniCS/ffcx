# Copyright (C) 2009-2017 Anders Logg and Martin Sandve Alnæs
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

# Note: Much of the code in this file is a direct translation
# from the old implementation in FFC, although some improvements
# have been made to the generated code.

import logging

import ffcx.codegeneration.finite_element_template as ufcx_finite_element
import ufl

logger = logging.getLogger("ffcx")
index_type = "int"


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
    d["element_type"] = ir.element_type
    d["space_dimension"] = ir.space_dimension
    d["value_rank"] = len(ir.value_shape)
    d["value_size"] = ufl.product(ir.value_shape)
    d["reference_value_rank"] = len(ir.reference_value_shape)
    d["reference_value_size"] = ufl.product(ir.reference_value_shape)
    d["degree"] = ir.degree
    d["family"] = f"\"{ir.family}\""
    d["num_sub_elements"] = ir.num_sub_elements
    d["block_size"] = ir.block_size
    d["discontinuous"] = "true" if ir.discontinuous else "false"

    if ir.lagrange_variant is None:
        d["lagrange_variant"] = -1
    else:
        d["lagrange_variant"] = int(ir.lagrange_variant)

    if ir.basix_family is None:
        d["basix_family"] = -1
    else:
        d["basix_family"] = int(ir.basix_family)
    if ir.basix_cell is None:
        d["basix_cell"] = -1
    else:
        d["basix_cell"] = int(ir.basix_cell)

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

    if len(ir.sub_elements) > 0:
        d["sub_elements"] = f"sub_elements_{ir.name}"
        d["sub_elements_init"] = L.ArrayDecl(
            "ufcx_finite_element*", f"sub_elements_{ir.name}",
            values=[L.AddressOf(L.Symbol(el)) for el in ir.sub_elements], sizes=len(ir.sub_elements))
    else:
        d["sub_elements"] = "NULL"
        d["sub_elements_init"] = ""

    # Check that no keys are redundant or have been missed
    from string import Formatter
    fieldnames = [
        fname for _, fname, _, _ in Formatter().parse(ufcx_finite_element.factory) if fname
    ]
    assert set(fieldnames) == set(
        d.keys()), "Mismatch between keys in template and in formattting dict"

    # Format implementation code
    implementation = ufcx_finite_element.factory.format_map(d)

    # Format declaration
    declaration = ufcx_finite_element.declaration.format(factory_name=ir.name)

    return declaration, implementation
