# Copyright (C) 2015-2020 Martin Sandve Alnæs and Chris Richardson
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import logging

import ffcx.codegeneration.coordinate_mapping_template as ufc_coordinate_mapping
import ffcx.codegeneration.C.cnodes as L

logger = logging.getLogger("ffcx")


def generator(ir, parameters):
    """Generate UFC code for a coordinate mapping."""

    logger.info("Generating code for coordinate mapping:")
    logger.info(f"--- cell shape: {ir.cell_shape}")
    logger.info(f"--- gdim: {ir.geometric_dimension}")
    logger.info(f"--- tdim: {ir.topological_dimension}")
    logger.info(f"--- name: {ir.name}")
    logger.info(f"--- scalar dofmap name: {ir.scalar_dofmap}")

    d = {}

    # Attributes
    d["factory_name"] = ir.name
    d["prefix"] = ir.prefix
    d["signature"] = f"\"{ir.signature}\""
    d["geometric_dimension"] = ir.geometric_dimension
    d["topological_dimension"] = ir.topological_dimension
    d["is_affine"] = 1 if ir.is_affine else 0
    d["cell_shape"] = ir.cell_shape
    d["scalar_dofmap"] = L.AddressOf(L.Symbol(ir.scalar_dofmap))

    d["family"] = f"\"{ir.coordinate_element_family}\""
    d["degree"] = ir.coordinate_element_degree

    # Check that no keys are redundant or have been missed
    from string import Formatter
    fields = [
        fname for _, fname, _, _ in Formatter().parse(ufc_coordinate_mapping.factory) if fname
    ]
    assert set(fields) == set(
        d.keys()), "Mismatch between keys in template and in formattting dict."

    # Format implementation code
    implementation = ufc_coordinate_mapping.factory.format_map(d)

    # Format declaration
    declaration = ufc_coordinate_mapping.declaration.format(factory_name=ir.name, prefix=ir.prefix)

    return declaration, implementation
