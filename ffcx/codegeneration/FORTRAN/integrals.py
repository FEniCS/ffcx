# Copyright (C) 2015-2021 Martin Sandve Alnæs, Michal Habera, Igor Baratta
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import logging

from ffcx.codegeneration.integral_generator import IntegralGenerator
from ffcx.codegeneration.FORTRAN import integrals_template as ufcx_integrals
from ffcx.codegeneration.backend import FFCXBackend
from ffcx.codegeneration.utils import cdtype_to_numpy
from ffcx.codegeneration.FORTRAN.fortran_implementation import FortranFormatter


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
    FF = FortranFormatter(options["scalar_type"])
    body = FF.c_format(parts)
    decl = FF.declarations
    body = decl + body

    # Generate generic FFCx code snippets and add specific parts
    code = {}
    code["class_type"] = ir.integral_type + "_integral"
    code["name"] = ir.name

    vals = ", ".join("1" if i else "0" for i in ir.enabled_coefficients)
    code["enabled_coefficients"] = f"(/ {vals} /)"

    code["additional_includes_set"] = set()  # FIXME: Get this out of code[]
    code["tabulate_tensor"] = body

    if options["tabulate_tensor_void"]:
        code["tabulate_tensor"] = ""

    scalar_types = {
        "double": "DOUBLE PRECISION",
        "float": "REAL",
        "_Complex double": "COMPLEX",
    }
    scalar_type = scalar_types.get(options["scalar_type"])
    value_type = scalar_type

    implementation = ufcx_integrals.factory.format(
        factory_name=factory_name,
        enabled_coefficients=code["enabled_coefficients"],
        tabulate_tensor=code["tabulate_tensor"],
        needs_facet_permutations="true" if ir.needs_facet_permutations else "false",
        scalar_type=scalar_type,
        geom_type=value_type,
        np_scalar_type=cdtype_to_numpy(options["scalar_type"]),
        coordinate_element=ir.coordinate_element,
    )

    return declaration, implementation
