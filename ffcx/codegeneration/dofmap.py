# Copyright (C) 2009-2018 Anders Logg, Martin Sandve Alnæs and Garth N. Wells
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# Note: Most of the code in this file is a direct translation from the
# old implementation in FFC

import logging

import ffcx.codegeneration.dofmap_template as ufcx_dofmap

logger = logging.getLogger("ffcx")


def tabulate_entity_dofs(L, entity_dofs, num_dofs_per_entity):
    # Output argument array
    dofs = L.Symbol("dofs")

    # Input arguments
    d = L.Symbol("d")
    i = L.Symbol("i")

    # TODO: Removed check for (d <= tdim + 1)
    tdim = len(num_dofs_per_entity) - 1

    # Generate cases for each dimension:
    all_cases = []
    for dim in range(tdim + 1):

        # Ignore if no entities for this dimension
        if num_dofs_per_entity[dim] == 0:
            continue

        # Generate cases for each mesh entity
        cases = []
        for entity in range(len(entity_dofs[dim])):
            casebody = []
            for (j, dof) in enumerate(entity_dofs[dim][entity]):
                casebody += [L.Assign(dofs[j], dof)]
            cases.append((entity, L.StatementList(casebody)))

        # Generate inner switch
        # TODO: Removed check for (i <= num_entities-1)
        inner_switch = L.Switch(i, cases, autoscope=False)
        all_cases.append((dim, inner_switch))

    if all_cases:
        return L.Switch(d, all_cases, autoscope=False)
    else:
        return L.NoOp()


def generator(ir, parameters):
    """Generate UFC code for a dofmap."""
    logger.info("Generating code for dofmap:")
    logger.info(f"--- num element support dofs: {ir.num_element_support_dofs}")
    logger.info(f"--- name: {ir.name}")

    d = {}

    # Attributes
    d["factory_name"] = ir.name
    d["signature"] = f"\"{ir.signature}\""
    d["num_global_support_dofs"] = ir.num_global_support_dofs
    d["num_element_support_dofs"] = ir.num_element_support_dofs
    d["num_sub_dofmaps"] = ir.num_sub_dofmaps

    import ffcx.codegeneration.C.cnodes as L

    num_entity_dofs = ir.num_entity_dofs + [0, 0, 0, 0]
    num_entity_dofs = num_entity_dofs[:4]
    d["num_entity_dofs"] = f"num_entity_dofs_{ir.name}"
    d["num_entity_dofs_init"] = L.ArrayDecl("int", f"num_entity_dofs_{ir.name}",
                                            values=num_entity_dofs, sizes=4)

    num_entity_closure_dofs = ir.num_entity_closure_dofs + [0, 0, 0, 0]
    num_entity_closure_dofs = num_entity_closure_dofs[:4]
    d["num_entity_closure_dofs"] = f"num_entity_closure_dofs_{ir.name}"
    d["num_entity_closure_dofs_init"] = L.ArrayDecl("int", f"num_entity_closure_dofs_{ir.name}",
                                                    values=num_entity_closure_dofs, sizes=4)

    d["block_size"] = ir.block_size

    # Functions
    d["tabulate_entity_dofs"] = tabulate_entity_dofs(L, *ir.tabulate_entity_dofs)
    d["tabulate_entity_closure_dofs"] = tabulate_entity_dofs(L, *ir.tabulate_entity_closure_dofs)

    if len(ir.sub_dofmaps) > 0:
        d["sub_dofmaps_initialization"] = L.ArrayDecl(
            "ufcx_dofmap*", f"sub_dofmaps_{ir.name}",
            values=[L.AddressOf(L.Symbol(dofmap)) for dofmap in ir.sub_dofmaps], sizes=len(ir.sub_dofmaps))
        d["sub_dofmaps"] = f"sub_dofmaps_{ir.name}"
    else:
        d["sub_dofmaps_initialization"] = ""
        d["sub_dofmaps"] = "NULL"

    # Check that no keys are redundant or have been missed
    from string import Formatter
    fields = [fname for _, fname, _, _ in Formatter().parse(ufcx_dofmap.factory) if fname]
    # Remove square brackets from any field names
    fields = [f.split("[")[0] for f in fields]
    assert set(fields) == set(
        d.keys()), "Mismatch between keys in template and in formattting dict."

    # Format implementation code
    implementation = ufcx_dofmap.factory.format_map(d)

    # Format declaration
    declaration = ufcx_dofmap.declaration.format(factory_name=ir.name)

    return declaration, implementation
