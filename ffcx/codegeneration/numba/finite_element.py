# Copyright (C) 2009-2022 Anders Logg, Martin Sandve Alnæs, Matthew Scroggs
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

# Note: Much of the code in this file is a direct translation
# from the old implementation in FFC, although some improvements
# have been made to the generated code.

import logging

import ffcx.codegeneration.numba.basix_custom_element_template as ufcx_basix_custom_finite_element
import ffcx.codegeneration.numba.finite_element_template as ufcx_finite_element
import ufl

logger = logging.getLogger("ffcx")
index_type = "int"


def generator(ir, options):
    """Generate UFC code for a finite element."""
    logger.info("Generating code for finite element:")
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
    d["num_sub_elements"] = ir.num_sub_elements
    d["block_size"] = ir.block_size
    d["discontinuous"] = str(ir.discontinuous)

    if ir.lagrange_variant is None:
        d["lagrange_variant"] = -1
    else:
        d["lagrange_variant"] = int(ir.lagrange_variant)

    if ir.dpc_variant is None:
        d["dpc_variant"] = -1
    else:
        d["dpc_variant"] = int(ir.dpc_variant)

    if ir.basix_family is None:
        d["basix_family"] = -1
    else:
        d["basix_family"] = int(ir.basix_family)
    if ir.basix_cell is None:
        d["basix_cell"] = -1
    else:
        d["basix_cell"] = int(ir.basix_cell)

    d["value_shape"] = f"[{', '.join(str(i) for i in ir.value_shape)}]"
    d["reference_value_shape"] = f"[{', '.join(str(i) for i in ir.reference_value_shape)}]"

    sub_elements = ", ".join(el for el in ir.sub_elements)
    d["sub_elements"] = f"[{sub_elements}]"

    d["custom_element"] = "None"
    if ir.custom_element is not None:
        d["custom_element"] = generate_custom_element(f"custom_element_{ir.name}", ir.custom_element)

    # Check that no keys are redundant or have been missed
    from string import Formatter
    fieldnames = [
        fname for _, fname, _, _ in Formatter().parse(ufcx_finite_element.factory) if fname
    ]
    assert set(fieldnames) == set(
        d.keys()), "Mismatch between keys in template and in formatting dict"

    # Format implementation code
    implementation = ufcx_finite_element.factory.format_map(d)

    # Format declaration
    declaration = ufcx_finite_element.declaration.format(factory_name=ir.name)

    return declaration, implementation


def generate_custom_element(name, ir):
    d = {}
    d["factory_name"] = name
    d["cell_type"] = int(ir.cell_type)
    d["polyset_type"] = int(ir.polyset_type)
    d["map_type"] = int(ir.map_type)
    d["sobolev_space"] = int(ir.sobolev_space)
    d["highest_complete_degree"] = ir.highest_complete_degree
    d["highest_degree"] = ir.highest_degree
    d["discontinuous"] = "true" if ir.discontinuous else "false"
    d["interpolation_nderivs"] = ir.interpolation_nderivs
    d["value_shape_length"] = len(ir.value_shape)
    if len(ir.value_shape) > 0:
        d["value_shape"] = f"value_shape_{name}"
        values = ', '.join(str(i) for i in ir.value_shape)
        d["value_shape_init"] = f"int value_shape_{name}[] = {{{values}}};"
    else:
        d["value_shape"] = "NULL"
        d["value_shape_init"] = ""

    d["wcoeffs_rows"] = ir.wcoeffs.shape[0]
    d["wcoeffs_cols"] = ir.wcoeffs.shape[1]
    d["wcoeffs"] = f"wcoeffs_{name}"
    d["wcoeffs_init"] = f"double wcoeffs_{name}[{ir.wcoeffs.shape[0] * ir.wcoeffs.shape[1]}] = "
    d["wcoeffs_init"] += "{" + ", ".join([f" {i}" for row in ir.wcoeffs for i in row]) + "};"

    npts = []
    x = []
    for entity in ir.x:
        for points in entity:
            npts.append(points.shape[0])
            for row in points:
                for i in row:
                    x.append(i)
    d["npts"] = f"npts_{name}"
    d["npts_init"] = f"int npts_{name}[{len(npts)}] = "
    d["npts_init"] += "{" + ",".join([f" {i}" for i in npts]) + "};"
    d["x"] = f"x_{name}"
    d["x_init"] = f"double x_{name}[{len(x)}] = "
    d["x_init"] += "{" + ",".join([f" {i}" for i in x]) + "};"
    ndofs = []
    M = []
    for entity in ir.M:
        for mat4d in entity:
            ndofs.append(mat4d.shape[0])
            for mat3d in mat4d:
                for mat2d in mat3d:
                    for row in mat2d:
                        for i in row:
                            M.append(i)

    d["ndofs"] = f"ndofs_{name}"
    d["ndofs_init"] = f"int ndofs_{name}[{len(ndofs)}] = "
    d["ndofs_init"] += "{" + ",".join([f" {i}" for i in ndofs]) + "};"
    d["M"] = f"M_{name}"
    d["M_init"] = f"double M_{name}[{len(M)}] = "
    d["M_init"] += "{" + ",".join([f" {i}" for i in M]) + "};"

    # Check that no keys are redundant or have been missed
    from string import Formatter
    fieldnames = [
        fname for _, fname, _, _ in Formatter().parse(ufcx_basix_custom_finite_element.factory) if fname
    ]
    assert set(fieldnames) == set(
        d.keys()), "Mismatch between keys in template and in formatting dict"

    # Format implementation code
    implementation = ufcx_basix_custom_finite_element.factory.format_map(d)

    return implementation
