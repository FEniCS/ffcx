# Copyright (C) 2009-2025 Anders Logg, Martin Sandve Alnæs, Chris Richardson and Paul T. Kühner
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

# Note: Most of the code in this file is a direct translation from the
# old implementation in FFC
"""Template for form output."""

import logging

import numpy as np
from numpy import typing as npt

from ffcx.codegeneration.common import template_keys
from ffcx.codegeneration.numba import form_template

logger = logging.getLogger("ffcx")


def generator(ir, options: dict[str, int | float | npt.DTypeLike]) -> tuple[str, str]:
    """Generate UFC code for a form."""
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

        d["original_coefficient_position_init"] = (
            f"original_coefficient_position_{ir.name} = [{values}]"
        )
        d["original_coefficient_positions"] = f"original_coefficient_position_{ir.name}"
    else:
        d["original_coefficient_position_init"] = ""
        d["original_coefficient_positions"] = "None"

    if len(ir.coefficient_names) > 0:
        values = ", ".join(f'"{name}"' for name in ir.coefficient_names)
        d["coefficient_names_init"] = f"coefficient_names_{ir.name} = [{values}]"
        d["coefficient_names"] = f"coefficient_names_{ir.name}"
    else:
        d["coefficient_names_init"] = ""
        d["coefficient_names"] = "None"

    d["num_constants"] = ir.num_constants
    if ir.num_constants > 0:
        d["constant_ranks_init"] = f"constant_ranks_{ir.name} = [{str(ir.constant_ranks)[1:-1]}]"
        d["constant_ranks"] = f"constant_ranks_{ir.name}"

        shapes = [
            f"constant_shapes_{ir.name}_{i} = [{str(shape)[1:-1]}]"
            for i, shape in enumerate(ir.constant_shapes)
            if len(shape) > 0
        ]
        names = [f"constant_shapes_{ir.name}_{i}" for i in range(ir.num_constants)]
        shapes1 = f"constant_shapes_{ir.name} = ["
        for rank, name in zip(ir.constant_ranks, names):
            if rank > 0:
                shapes1 += f"{name},\n"
            else:
                shapes1 += "None,\n"
        shapes1 += "]"
        shapes.append(shapes1)

        d["constant_shapes_init"] = "\n".join(shapes)
        d["constant_shapes"] = f"constant_shapes_{ir.name}"
    else:
        d["constant_ranks_init"] = ""
        d["constant_ranks"] = "None"
        d["constant_shapes_init"] = ""
        d["constant_shapes"] = "None"

    if len(ir.constant_names) > 0:
        values = ", ".join(f'"{name}"' for name in ir.constant_names)
        d["constant_names_init"] = f"constant_names_{ir.name} = [{values}]"
        d["constant_names"] = f"constant_names_{ir.name}"
    else:
        d["constant_names_init"] = ""
        d["constant_names"] = "None"

    if len(ir.finite_element_hashes) > 0:
        d["finite_element_hashes"] = f"finite_element_hashes_{ir.name}"
        values = ", ".join(f"{0 if el is None else el}" for el in ir.finite_element_hashes)
        d["finite_element_hashes_init"] = f"finite_element_hashes_{ir.name} = [{values}]"
    else:
        d["finite_element_hashes"] = "None"
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
            unsorted_integrals += [name]
            unsorted_ids += [id]
            unsorted_domains += [domains]

        id_sort = np.argsort(unsorted_ids)
        integrals += [unsorted_integrals[i] for i in id_sort]
        integral_ids += [unsorted_ids[i] for i in id_sort]
        integral_domains += [unsorted_domains[i] for i in id_sort]

        integral_offsets.append(sum(len(d) for d in integral_domains))

    if len(integrals) > 0:
        values = ", ".join(
            [
                f"{i}_{domain.name}"
                for i, domains in zip(integrals, integral_domains)
                for domain in domains
            ]
        )
        d["form_integrals_init"] = f"form_integrals_{ir.name} = [{values}]"
        d["form_integrals"] = f"form_integrals_{ir.name}"
        values = ", ".join(
            f"{i}" for i, domains in zip(integral_ids, integral_domains) for _ in domains
        )
        d["form_integral_ids_init"] = f"form_integral_ids_{ir.name} = [{values}]"
        d["form_integral_ids"] = f"form_integral_ids_{ir.name}"
    else:
        d["form_integrals_init"] = ""
        d["form_integrals"] = "None"
        d["form_integral_ids_init"] = ""
        d["form_integral_ids"] = "None"

    values = ", ".join(str(i) for i in integral_offsets)
    d["form_integral_offsets_init"] = f"form_integral_offsets_{ir.name} = [{values}]"

    # Format implementation code
    assert set(d.keys()) == template_keys(form_template.factory)
    implementation = form_template.factory.format_map(d)

    return "", implementation
