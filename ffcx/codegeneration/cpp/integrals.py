# Copyright (C) 2015-2021 Martin Sandve Alnæs, Michal Habera, Igor Baratta
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Integral generation."""

import logging

from ffcx.codegeneration.backend import FFCXBackend
from ffcx.codegeneration.cpp import integrals_template as ufcx_integrals
from ffcx.codegeneration.cpp.cpp_implementation import CppFormatter
from ffcx.codegeneration.integral_generator import IntegralGenerator

logger = logging.getLogger("ffcx")


def generator(ir, options):
    """Generate code for an integral."""
    logger.info("Generating code for integral:")
    logger.info(f"--- type: {ir.integral_type}")
    logger.info(f"--- name: {ir.name}")

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
    CF = CppFormatter(options["scalar_type"])
    body = CF.c_format(parts)

    # Generate generic FFCx code snippets and add specific parts
    code = {}
    code["class_type"] = ir.integral_type + "_integral"
    code["name"] = ir.name

    vals = ", ".join("true" if i else "false" for i in ir.enabled_coefficients)
    code["enabled_coefficients"] = f"{{{vals}}}"

    code["additional_includes_set"] = set()  # FIXME: Get this out of code[]
    code["tabulate_tensor"] = body

    implementation = ufcx_integrals.factory.format(
        factory_name=factory_name,
        enabled_coefficients=code["enabled_coefficients"],
        tabulate_tensor=code["tabulate_tensor"],
        needs_facet_permutations="true" if ir.needs_facet_permutations else "false",
        scalar_type=options["scalar_type"],
        geom_type=options["scalar_type"],
        np_scalar_type=options["scalar_type"],
        coordinate_element=ir.coordinate_element,
    )

    return declaration + implementation, ""
