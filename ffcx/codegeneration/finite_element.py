# Copyright (C) 2009-2017 Anders Logg and Martin Sandve Alnæs
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

# Note: Much of the code in this file is a direct translation
# from the old implementation in FFC, although some improvements
# have been made to the generated code.

import logging

import ffcx.codegeneration.finite_element_template as ufc_finite_element
import ufl

logger = logging.getLogger("ffcx")
index_type = "int"


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

    if ir.tabulate_entity_dofs is None:
        d["num_entity_dofs"] = [-1, -1, -1, -1]
        d["tabulate_entity_dofs"] = L.NoOp()
        d["num_entity_closure_dofs"] = [-1, -1, -1, -1]
        d["tabulate_entity_closure_dofs"] = L.NoOp()
    else:
        d["num_entity_dofs"] = ir.num_entity_dofs + [0, 0, 0, 0]
        d["tabulate_entity_dofs"] = tabulate_entity_dofs(L, *ir.tabulate_entity_dofs)
        d["num_entity_closure_dofs"] = ir.num_entity_closure_dofs + [0, 0, 0, 0]
        d["tabulate_entity_closure_dofs"] = tabulate_entity_dofs(L, *ir.tabulate_entity_closure_dofs)

    if len(ir.sub_elements) > 0:
        d["sub_elements"] = f"sub_elements_{ir.name}"
        d["sub_elements_init"] = L.ArrayDecl(
            "ufc_finite_element*", f"sub_elements_{ir.name}",
            values=[L.AddressOf(L.Symbol(el)) for el in ir.sub_elements], sizes=len(ir.sub_elements))
    else:
        d["sub_elements"] = "NULL"
        d["sub_elements_init"] = ""

    # Check that no keys are redundant or have been missed
    from string import Formatter
    fieldnames = [
        fname.split("[")[0] for _, fname, _, _ in Formatter().parse(ufc_finite_element.factory) if fname
    ]

    assert set(fieldnames) == set(
        d.keys()), "Mismatch between keys in template and in formattting dict"

    # Format implementation code
    implementation = ufc_finite_element.factory.format_map(d)

    # Format declaration
    declaration = ufc_finite_element.declaration.format(factory_name=ir.name)

    return declaration, implementation
