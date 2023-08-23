# Copyright (C) 2009-2018 Anders Logg, Martin Sandve Alnæs and Garth N. Wells
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# Note: Most of the code in this file is a direct translation from the
# old implementation in FFC

import logging

import ffcx.codegeneration.C.dofmap_template as ufcx_dofmap

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

    import ffcx.codegeneration.C.cnodes as L

    flattened_entity_dofs = []
    entity_dof_offsets = [0]
    for dim in ir.entity_dofs:
        for ent in dim:
            for v in ent:
                flattened_entity_dofs.append(v)
            entity_dof_offsets.append(len(flattened_entity_dofs))
    d["entity_dofs"] = f"entity_dofs_{ir.name}"
    d["entity_dofs_init"] = L.ArrayDecl(
        "int",
        f"entity_dofs_{ir.name}",
        values=flattened_entity_dofs,
        sizes=len(flattened_entity_dofs),
    )
    d["entity_dof_offsets"] = f"entity_dof_offsets_{ir.name}"
    d["entity_dof_offsets_init"] = L.ArrayDecl(
        "int",
        f"entity_dof_offsets_{ir.name}",
        values=entity_dof_offsets,
        sizes=len(entity_dof_offsets),
    )

    # Closure
    flattened_entity_closure_dofs = []
    entity_closure_dof_offsets = [0]
    for dim in ir.entity_closure_dofs:
        for ent in dim:
            for v in ent:
                flattened_entity_closure_dofs.append(v)
            entity_closure_dof_offsets.append(len(flattened_entity_closure_dofs))
    d["entity_closure_dofs"] = f"entity_closure_dofs_{ir.name}"
    d["entity_closure_dofs_init"] = L.ArrayDecl(
        "int",
        f"entity_closure_dofs_{ir.name}",
        values=flattened_entity_closure_dofs,
        sizes=len(flattened_entity_closure_dofs),
    )
    d["entity_closure_dof_offsets"] = f"entity_closure_dof_offsets_{ir.name}"
    d["entity_closure_dof_offsets_init"] = L.ArrayDecl(
        "int",
        f"entity_closure_dof_offsets_{ir.name}",
        values=entity_closure_dof_offsets,
        sizes=len(entity_closure_dof_offsets),
    )

    d["block_size"] = ir.block_size

    if len(ir.sub_dofmaps) > 0:
        d["sub_dofmaps_initialization"] = L.ArrayDecl(
            "ufcx_dofmap*",
            f"sub_dofmaps_{ir.name}",
            values=[L.AddressOf(L.Symbol(dofmap)) for dofmap in ir.sub_dofmaps],
            sizes=len(ir.sub_dofmaps),
        )
        d["sub_dofmaps"] = f"sub_dofmaps_{ir.name}"
    else:
        d["sub_dofmaps_initialization"] = ""
        d["sub_dofmaps"] = "NULL"

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

    # Format declaration
    declaration = ufcx_dofmap.declaration.format(factory_name=ir.name)

    return declaration, implementation
