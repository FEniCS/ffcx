# Copyright (C) 2015-2021 Martin Sandve Alnæs, Michal Habera, Igor Baratta
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import logging

from ffcx.codegeneration.integral_generator import IntegralGenerator
from ffcx.codegeneration.C import integrals_template as ufcx_integrals
from ffcx.codegeneration.backend import FFCXBackend
from ffcx.codegeneration.C.format_lines import format_indented_lines
from ffcx.naming import cdtype_to_numpy, scalar_to_value_type

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
    body = format_indented_lines(parts.cs_format(ir.precision), 1)

    # Generate generic FFCx code snippets and add specific parts
    code = {}
    code["class_type"] = ir.integral_type + "_integral"
    code["name"] = ir.name
    code["members"] = ""
    code["constructor"] = ""
    code["constructor_arguments"] = ""
    code["initializer_list"] = ""
    code["destructor"] = ""

    L = backend.language
    if len(ir.enabled_coefficients) > 0:
        code["enabled_coefficients_init"] = L.ArrayDecl(
            "bool", f"enabled_coefficients_{ir.name}", values=ir.enabled_coefficients,
            sizes=len(ir.enabled_coefficients))
        code["enabled_coefficients"] = f"enabled_coefficients_{ir.name}"
    else:
        code["enabled_coefficients_init"] = ""
        code["enabled_coefficients"] = L.Null()

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
        coordinate_element=L.AddressOf(L.Symbol(ir.coordinate_element)))

    return declaration, implementation
