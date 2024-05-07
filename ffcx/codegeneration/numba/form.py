# Copyright (C) 2009-2017 Anders Logg and Martin Sandve Alnæs
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

# Note: Most of the code in this file is a direct translation from the
# old implementation in FFC
"""Template for form output."""

import logging

from ffcx.codegeneration.numba import form_template

logger = logging.getLogger("ffcx")


def generator(ir, options):
    """Generate UFC code for a form."""
    logger.info("Generating code for form:")
    logger.info(f"--- rank: {ir.rank}")
    logger.info(f"--- name: {ir.name}")

    d = {}
    d["factory_name"] = ir.name
    d["name_from_uflfile"] = ir.name_from_uflfile
    d["signature"] = f'"{ir.signature}"'
    d["rank"] = ir.rank
    d["num_coefficients"] = ir.num_coefficients
    d["num_constants"] = ir.num_constants

    orig_coeff = ', '.join(str(i) for i in ir.original_coefficient_position)
    d["original_coefficient_position"] = f"[{orig_coeff}]"

    cnames = ir.coefficient_names
    assert ir.num_coefficients == len(cnames)
    names = ", ".join(f'"{name}"' for name in cnames)
    d["coefficient_name_map"] = f"[{names}]"

    cstnames = ir.constant_names
    names = ", ".join(f'"{name}"' for name in cstnames)
    d["constant_name_map"] = f"[{names}]"

    integrals = []
    integral_ids = []
    integral_offsets = [0]
    for itg_type in ("cell", "exterior_facet", "interior_facet"):
        integrals += ir.integral_names[itg_type]
        integral_ids += ir.subdomain_ids[itg_type]
        integral_offsets.append(len(integrals))

    integrals_str = ", ".join(itg for itg in integrals)
    integral_ids_str = ", ".join(str(i) for i in integral_ids)
    d["form_integrals"] = f"[{integrals_str}]"
    d["form_integral_ids"] = f"[{integral_ids_str}]"
    offsets = ", ".join(str(i) for i in integral_offsets)
    d["form_integral_offsets"] = f"[{offsets}]"

    # Check that no keys are redundant or have been missed
    from string import Formatter

    fields = [
        fname for _, fname, _, _ in Formatter().parse(form_template.factory) if fname
    ]

    for f in fields:
        if f not in d.keys():
            print(f, "not in d.keys()")

    for f in d.keys():
        if f not in fields:
            print(f, "not in fields")

    if set(fields) != set(d.keys()):
        print("Mismatch between keys in template and in formatting dict")

    # Format implementation code
    implementation = form_template.factory.format_map(d)

    return "", implementation
