# Copyright (C) 2009-2018 Anders Logg, Martin Sandve Alnæs and Garth N. Wells
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# Note: Most of the code in this file is a direct translation from the
# old implementation in FFC

import logging

import ffcx.codegeneration.FORTRAN.dofmap_template as ufcx_dofmap

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

    flattened_entity_dofs = []
    entity_dof_offsets = [0]
    for dim in ir.entity_dofs:
        for ent in dim:
            for v in ent:
                flattened_entity_dofs.append(v)
            entity_dof_offsets.append(len(flattened_entity_dofs))
    values = ", ".join(str(i) for i in flattened_entity_dofs)
    d["entity_dofs"] = f"(/{values}/)"
    values = ", ".join(str(i) for i in entity_dof_offsets)
    d["entity_dof_offsets"] = f"(/{values}/)"

    # Closure
    flattened_entity_closure_dofs = []
    entity_closure_dof_offsets = [0]
    for dim in ir.entity_closure_dofs:
        for ent in dim:
            for v in ent:
                flattened_entity_closure_dofs.append(v)
            entity_closure_dof_offsets.append(len(flattened_entity_closure_dofs))
    values = ", ".join(str(i) for i in flattened_entity_closure_dofs)
    d["entity_closure_dofs"] = f"(/{values}/)"
    values = ", ".join(str(i) for i in entity_closure_dof_offsets)
    d["entity_closure_dof_offsets"] = f"(/{values}/)"

    d["block_size"] = ir.block_size

    if len(ir.sub_dofmaps) > 0:
        values = ", ".join(f"{dofmap}" for dofmap in ir.sub_dofmaps)
        d["sub_dofmaps"] = f"(/{values}/)"
    else:
        d["sub_dofmaps"] = "NULL"

    # Check that no keys are redundant or have been missed
    # from string import Formatter

    # fields = [
    #     fname for _, fname, _, _ in Formatter().parse(ufcx_dofmap.factory) if fname
    # ]
    # # Remove square brackets from any field names
    # fields = [f.split("[")[0] for f in fields]
    # assert set(fields) == set(
    #     d.keys()
    # ), "Mismatch between keys in template and in formatting dict."

    # Format implementation code
    implementation = ufcx_dofmap.factory.format_map(d)

    # Format declaration
    declaration = ufcx_dofmap.declaration.format(factory_name=ir.name)

    return declaration, implementation
