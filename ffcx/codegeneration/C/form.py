# Copyright (C) 2009-2017 Anders Logg and Martin Sandve Alnæs
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# Modified by Chris Richardson and Jørgen S. Dokken 2023
#
# Note: Most of the code in this file is a direct translation from the
# old implementation in FFC
"""Generate UFC code for a form."""

from __future__ import annotations

import logging

import numpy as np

from ffcx.codegeneration.C import form_template
from ffcx.ir.representation import FormIR

logger = logging.getLogger("ffcx")


def generator(ir: FormIR, options):
    """Generate UFCx code for a form."""
    logger.info("Generating code for form:")
    logger.info(f"--- rank: {ir.rank}")
    logger.info(f"--- name: {ir.name}")

    d: dict[str, int | str] = {}
    d["factory_name"] = ir.name
    d["name_from_uflfile"] = ir.name_from_uflfile
    d["signature"] = f'"{ir.signature}"'
    d["rank"] = ir.rank
    d["num_coefficients"] = ir.num_coefficients

    if len(ir.original_coefficient_positions) > 0:
        values = ", ".join(str(i) for i in ir.original_coefficient_positions)
        sizes = len(ir.original_coefficient_positions)

        d["original_coefficient_position_init"] = (
            f"int original_coefficient_position_{ir.name}[{sizes}] = {{{values}}};"
        )
        d["original_coefficient_positions"] = f"original_coefficient_position_{ir.name}"
    else:
        d["original_coefficient_position_init"] = ""
        d["original_coefficient_positions"] = "NULL"

    if len(ir.coefficient_names) > 0:
        values = ", ".join(f'"{name}"' for name in ir.coefficient_names)
        sizes = len(ir.coefficient_names)
        d["coefficient_names_init"] = (
            f"static const char* coefficient_names_{ir.name}[{sizes}] = {{{values}}};"
        )
        d["coefficient_names"] = f"coefficient_names_{ir.name}"
    else:
        d["coefficient_names_init"] = ""
        d["coefficient_names"] = "NULL"

    d["num_constants"] = ir.num_constants
    if ir.num_constants > 0:
        d["constant_ranks_init"] = (
            f"static const int constant_ranks_{ir.name}[{ir.num_constants}] = "
            f"{{{str(ir.constant_ranks)[1:-1]}}};"
        )
        d["constant_ranks"] = f"constant_ranks_{ir.name}"

        shapes = [
            f"static const int constant_shapes_{ir.name}_{i}[{len(shape)}] = "
            f"{{{str(shape)[1:-1]}}};"
            for i, shape in enumerate(ir.constant_shapes)
            if len(shape) > 0
        ]
        names = [f"constant_shapes_{ir.name}_{i}" for i in range(ir.num_constants)]
        shapes1 = f"static const int* constant_shapes_{ir.name}[{ir.num_constants}] = " + "{"
        for rank, name in zip(ir.constant_ranks, names):
            if rank > 0:
                shapes1 += f"{name},\n"
            else:
                shapes1 += "NULL,\n"
        shapes1 += "};"
        shapes.append(shapes1)

        d["constant_shapes_init"] = "\n".join(shapes)
        d["constant_shapes"] = f"constant_shapes_{ir.name}"
    else:
        d["constant_ranks_init"] = ""
        d["constant_ranks"] = "NULL"
        d["constant_shapes_init"] = ""
        d["constant_shapes"] = "NULL"

    if len(ir.constant_names) > 0:
        values = ", ".join(f'"{name}"' for name in ir.constant_names)
        sizes = len(ir.constant_names)
        d["constant_names_init"] = (
            f"static const char* constant_names_{ir.name}[{sizes}] = {{{values}}};"
        )
        d["constant_names"] = f"constant_names_{ir.name}"
    else:
        d["constant_names_init"] = ""
        d["constant_names"] = "NULL"

    if len(ir.finite_element_hashes) > 0:
        d["finite_element_hashes"] = f"finite_element_hashes_{ir.name}"
        values = ", ".join(
            f"UINT64_C({0 if el is None else el})" for el in ir.finite_element_hashes
        )
        sizes = len(ir.finite_element_hashes)
        d["finite_element_hashes_init"] = (
            f"uint64_t finite_element_hashes_{ir.name}[{sizes}] = {{{values}}};"
        )
    else:
        d["finite_element_hashes"] = "NULL"
        d["finite_element_hashes_init"] = ""

    integrals = []
    integral_ids = []
    integral_offsets = [0]
    integral_domains = []
    # Note: the order of this list is defined by the enum ufcx_integral_type in ufcx.h
    for itg_type in ("cell", "exterior_facet", "interior_facet", "vertex", "ridge"):
        unsorted_integrals = []
        unsorted_ids = []
        unsorted_domains = []
        for name, domains, id in zip(
            ir.integral_names[itg_type],
            ir.integral_domains[itg_type],
            ir.subdomain_ids[itg_type],
        ):
            unsorted_integrals += [f"&{name}"]
            unsorted_ids += [id]
            unsorted_domains += [domains]

        id_sort = np.argsort(unsorted_ids)
        integrals += [unsorted_integrals[i] for i in id_sort]
        integral_ids += [unsorted_ids[i] for i in id_sort]
        integral_domains += [unsorted_domains[i] for i in id_sort]

        integral_offsets.append(sum(len(d) for d in integral_domains))

    if len(integrals) > 0:
        sizes = sum(len(domains) for domains in integral_domains)
        values = ", ".join(
            [
                f"{i}_{domain.name}"
                for i, domains in zip(integrals, integral_domains)
                for domain in domains
            ]
        )
        d["form_integrals_init"] = (
            f"static ufcx_integral* form_integrals_{ir.name}[{sizes}] = {{{values}}};"
        )
        d["form_integrals"] = f"form_integrals_{ir.name}"
        values = ", ".join(
            f"{i}" for i, domains in zip(integral_ids, integral_domains) for _ in domains
        )
        d["form_integral_ids_init"] = f"int form_integral_ids_{ir.name}[{sizes}] = {{{values}}};"
        d["form_integral_ids"] = f"form_integral_ids_{ir.name}"
    else:
        d["form_integrals_init"] = ""
        d["form_integrals"] = "NULL"
        d["form_integral_ids_init"] = ""
        d["form_integral_ids"] = "NULL"

    sizes = len(integral_offsets)
    values = ", ".join(str(i) for i in integral_offsets)
    d["form_integral_offsets_init"] = (
        f"int form_integral_offsets_{ir.name}[{sizes}] = {{{values}}};"
    )

    # Check that no keys are redundant or have been missed
    from string import Formatter

    fields = [fname for _, fname, _, _ in Formatter().parse(form_template.factory) if fname]
    assert set(fields) == set(d.keys()), "Mismatch between keys in template and in formatting dict"

    # Format implementation code
    implementation = form_template.factory.format_map(d)

    # Format declaration
    declaration = form_template.declaration.format(
        factory_name=d["factory_name"], name_from_uflfile=d["name_from_uflfile"]
    )

    return declaration, implementation
