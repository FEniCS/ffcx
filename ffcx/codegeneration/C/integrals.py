# Copyright (C) 2015-2021 Martin Sandve Alnæs, Michal Habera, Igor Baratta
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import logging

from ffcx.codegeneration.integral_generator import IntegralGenerator
from ffcx.codegeneration.C import integrals_template as ufcx_integrals
from ffcx.codegeneration.backend import FFCXBackend
from ffcx.codegeneration.C.c_implementation import CFormatter
from ffcx.codegeneration.utils import cdtype_to_numpy, scalar_to_value_type

logger = logging.getLogger("ffcx")


def generator(ir, options):
    logger.info("Generating code for integral:")
    logger.info(f"--- type: {ir.integral_type}")
    logger.info(f"--- name: {ir.name}")

    """Generate code for an integral."""
    factory_name = ir.name

    # Format declaration
    declaration = ufcx_integrals.declaration.format(factory_name=factory_name)

    # Create FFCx C backend
    backend = FFCXBackend(ir, options)

    # Configure kernel generator
    ig = IntegralGenerator(ir, backend)

    # Generate code ast for the tabulate_tensor body
    parts = ig.generate()

    # Format code as string
    CF = CFormatter(options["scalar_type"])
    body = CF.c_format(parts)

    # Generate generic FFCx code snippets and add specific parts
    code = {}

    if len(ir.enabled_coefficients) > 0:
        values = ", ".join("1" if i else "0" for i in ir.enabled_coefficients)
        sizes = len(ir.enabled_coefficients)
        code["enabled_coefficients_init"] = f"bool enabled_coefficients_{ir.name}[{sizes}] = {{{values}}};"
        code["enabled_coefficients"] = f"enabled_coefficients_{ir.name}"
    else:
        code["enabled_coefficients_init"] = ""
        code["enabled_coefficients"] = "NULL"

    code["additional_includes_set"] = set()  # FIXME: Get this out of code[]
    code["tabulate_tensor"] = body

    implementation = ufcx_integrals.factory.format(
        factory_name=factory_name,
        enabled_coefficients=code["enabled_coefficients"],
        enabled_coefficients_init=code["enabled_coefficients_init"],
        tabulate_tensor=code["tabulate_tensor"],
        needs_facet_permutations="true" if ir.needs_facet_permutations else "false",
        scalar_type=options["scalar_type"],
        geom_type=scalar_to_value_type(options["scalar_type"]),
        np_scalar_type=cdtype_to_numpy(options["scalar_type"]),
        coordinate_element=f"&{ir.coordinate_element}")

    return declaration, implementation
