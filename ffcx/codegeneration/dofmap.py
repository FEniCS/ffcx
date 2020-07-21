# Copyright (C) 2009-2018 Anders Logg, Martin Sandve Alnæs and Garth N. Wells
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# Note: Most of the code in this file is a direct translation from the
# old implementation in FFC

import logging
import ffcx.codegeneration.dofmap_template as ufc_dofmap
from ffcx.codegeneration.utils import generate_return_new_switch

logger = logging.getLogger("ffcx")


def tabulate_entity_dofs(L, ir):
    entity_dofs, num_dofs_per_entity = ir.tabulate_entity_dofs

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


def create_sub_dofmap(L, ir):
    classnames = ir.create_sub_dofmap
    return generate_return_new_switch(L, "i", classnames)


def sub_dofmap_declaration(L, ir):
    classnames = set(ir.create_sub_dofmap)
    code = ""
    for name in classnames:
        code += "ufc_dofmap* create_{name}(void);\n".format(name=name)
    return code


def generator(ir, parameters):
    """Generate UFC code for a dofmap."""

    logger.info("Generating code for dofmap:")
    logger.info("--- num element support dofs: {}".format(ir.num_element_support_dofs))
    logger.info("--- name: {}".format(ir.name))

    d = {}

    # Attributes
    d["factory_name"] = ir.name
    d["signature"] = "\"{}\"".format(ir.signature)
    d["num_global_support_dofs"] = ir.num_global_support_dofs
    d["num_element_support_dofs"] = ir.num_element_support_dofs
    d["num_sub_dofmaps"] = ir.num_sub_dofmaps
    d["num_entity_dofs"] = ir.num_entity_dofs + [0, 0, 0, 0]
    d["block_size"] = ir.block_size

    num_perms = len(ir.base_permutations)
    if num_perms == 0:
        num_dofs = 0
    else:
        num_dofs = len(ir.base_permutations[0])

    bp = []
    for i, perm in enumerate(ir.base_permutations):
        for j, val in enumerate(perm):
            bp.append(str(val))
    d["base_permutations"] = ("static const int bp[" + str(num_perms * num_dofs) + "] = {"
                              + ",".join(bp) + "};\n  dofmap->base_permutations = bp;\n")
    d["size_base_permutations"] = num_perms * num_dofs

    import ffcx.codegeneration.C.cnodes as L

    # Functions
    d["tabulate_entity_dofs"] = tabulate_entity_dofs(L, ir)
    d["sub_dofmap_declaration"] = sub_dofmap_declaration(L, ir)
    d["create_sub_dofmap"] = create_sub_dofmap(L, ir)

    # Check that no keys are redundant or have been missed
    from string import Formatter
    fields = [fname for _, fname, _, _ in Formatter().parse(ufc_dofmap.factory) if fname]
    # Remove square brackets from any field names
    fields = [f.split("[")[0] for f in fields]
    assert set(fields) == set(
        d.keys()), "Mismatch between keys in template and in formattting dict."

    # Format implementation code
    implementation = ufc_dofmap.factory.format_map(d)

    # Format declaration
    declaration = ufc_dofmap.declaration.format(factory_name=ir.name)

    return declaration, implementation
