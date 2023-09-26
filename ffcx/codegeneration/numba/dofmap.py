# Copyright (C) 2009-2018 Anders Logg, Martin Sandve Alnæs and Garth N. Wells
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# Note: Most of the code in this file is a direct translation from the
# old implementation in FFC

import logging

import ffcx.codegeneration.numba.dofmap_template as ufcx_dofmap

logger = logging.getLogger("ffcx")


def generator(ir, options):
    """Generate UFC code for a dofmap."""
    logger.info("Generating code for dofmap:")
    logger.info(f"--- num element support dofs: {ir.num_element_support_dofs}")
    logger.info(f"--- name: {ir.name}")

    d = {}

    # Attributes
    d["factory_name"] = ir.name
    d["signature"] = f'"{ir.signature}"'
    d["num_global_support_dofs"] = ir.num_global_support_dofs
    d["num_element_support_dofs"] = ir.num_element_support_dofs
    d["num_sub_dofmaps"] = ir.num_sub_dofmaps
    d["entity_dofs"] = f"{ir.entity_dofs}"

    # Closure
    d["entity_closure_dofs"] = f"{ir.entity_closure_dofs}"

    d["block_size"] = ir.block_size

    values = ", ".join(dofmap for dofmap in ir.sub_dofmaps)
    d["sub_dofmaps"] = f"[{values}]"

    # Check that no keys are redundant or have been missed
    from string import Formatter

    fields = [
        fname for _, fname, _, _ in Formatter().parse(ufcx_dofmap.factory) if fname
    ]
    # Remove square brackets from any field names
    fields = [f.split("[")[0] for f in fields]
    assert set(fields) == set(
        d.keys()
    ), "Mismatch between keys in template and in formatting dict."

    # Format implementation code
    implementation = ufcx_dofmap.factory.format_map(d)
    return "", implementation
