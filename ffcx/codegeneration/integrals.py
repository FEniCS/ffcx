# Copyright (C) 2015-2021 Martin Sandve Alnæs, Michal Habera, Igor Baratta
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import logging

from ffcx.codegeneration.integral_generator import IntegralGenerator
from ffcx.codegeneration import integrals_template as ufcx_integrals
from ffcx.codegeneration.backend import FFCXBackend
from ffcx.codegeneration.C.format_lines import format_indented_lines
from ffcx.naming import cdtype_to_numpy, scalar_to_value_type
from ffcx.codegeneration.c_implementation import CFormatter


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
    code = CF.c_format(parts)
    body = format_indented_lines(code, 1)

    # Generate generic FFCx code snippets and add specific parts
    code = {}
    code["class_type"] = ir.integral_type + "_integral"
    code["name"] = ir.name

    if len(ir.enabled_coefficients) > 0:
        vals = ", ".join("1" if i else "0" for i in ir.enabled_coefficients)
        n = len(ir.enabled_coefficients)
        code[
            "enabled_coefficients_init"
        ] = f"bool enabled_coefficients_{ir.name}[{n}] = {{{vals}}};"

        code["enabled_coefficients"] = f"enabled_coefficients_{ir.name}"
    else:
        code["enabled_coefficients_init"] = ""
        code["enabled_coefficients"] = "NULL"

    code["additional_includes_set"] = set()  # FIXME: Get this out of code[]
    code["tabulate_tensor"] = body

    if options["tabulate_tensor_void"]:
        code["tabulate_tensor"] = ""

    implementation = ufcx_integrals.factory.format(
        factory_name=factory_name,
        enabled_coefficients=code["enabled_coefficients"],
        enabled_coefficients_init=code["enabled_coefficients_init"],
        tabulate_tensor=code["tabulate_tensor"],
        needs_facet_permutations="true" if ir.needs_facet_permutations else "false",
        scalar_type=options["scalar_type"],
        geom_type=scalar_to_value_type(options["scalar_type"]),
        np_scalar_type=cdtype_to_numpy(options["scalar_type"]),
        coordinate_element=ir.coordinate_element,
    )

    return declaration, implementation
