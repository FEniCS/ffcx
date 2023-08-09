# Copyright (C) 2009-2018 Anders Logg, Martin Sandve Alnæs and Garth N. Wells
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# Note: Most of the code in this file is a direct translation from the
# old implementation in FFC

import logging
import typing

import ffcx.codegeneration.C.dofmap_template as ufcx_dofmap

logger = logging.getLogger("ffcx")


def tabulate_entity_dofs(
    entity_dofs: typing.List[typing.List[typing.List[int]]],
    num_dofs_per_entity: typing.List[int],
):
    # TODO: Removed check for (d <= tdim + 1)
    tdim = len(num_dofs_per_entity) - 1

    # Generate cases for each dimension:
    cases = "switch (d)\n{\n"
    for dim in range(tdim + 1):
        # Ignore if no entities for this dimension
        if num_dofs_per_entity[dim] == 0:
            continue
        cases += f"  case {dim}:\n"
        cases += "   switch (i)\n   {\n"
        # Generate cases for each mesh entity
        for entity in range(len(entity_dofs[dim])):
            cases += f"  case {entity}:\n"
            for j, dof in enumerate(entity_dofs[dim][entity]):
                cases += f"     dofs[{j}] = {dof};\n"
            cases += "     break;\n"

        cases += "  }\n  break;\n"
    cases += "}\n"
    return cases


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

    num_entity_dofs = ir.num_entity_dofs + [0, 0, 0, 0]
    num_entity_dofs = num_entity_dofs[:4]
    d["num_entity_dofs"] = f"num_entity_dofs_{ir.name}"
    d[
        "num_entity_dofs_init"
    ] = f"int num_entity_dofs_{ir.name}[4] = {{{', '.join(str(i) for i in num_entity_dofs)}}};"

    num_entity_closure_dofs = ir.num_entity_closure_dofs + [0, 0, 0, 0]
    num_entity_closure_dofs = num_entity_closure_dofs[:4]
    d["num_entity_closure_dofs"] = f"num_entity_closure_dofs_{ir.name}"
    d[
        "num_entity_closure_dofs_init"
    ] = f"int num_entity_closure_dofs_{ir.name}[4] = {{{', '.join(str(i) for i in num_entity_closure_dofs)}}};"
    d["block_size"] = ir.block_size

    # Functions
    d["tabulate_entity_dofs"] = tabulate_entity_dofs(ir.entity_dofs, ir.num_entity_dofs)
    d["tabulate_entity_closure_dofs"] = tabulate_entity_dofs(
        ir.entity_closure_dofs, ir.num_entity_closure_dofs
    )

    if len(ir.sub_dofmaps) > 0:
        sub_dofmaps = ", ".join(f"&{dm}" for dm in ir.sub_dofmaps)
        d[
            "sub_dofmaps_initialization"
        ] = f"ufcx_dofmap* sub_dofmaps_{ir.name}[] = {{{sub_dofmaps}}};"
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
