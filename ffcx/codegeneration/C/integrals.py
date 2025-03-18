# Copyright (C) 2015-2021 Martin Sandve Alnæs, Michal Habera, Igor Baratta
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Generate UFC code for an integral."""

import logging
import sys

import basix
import numpy as np

from ffcx.codegeneration.backend import FFCXBackend
from ffcx.codegeneration.C import integrals_template as ufcx_integrals
from ffcx.codegeneration.C.c_implementation import CFormatter
from ffcx.codegeneration.integral_generator import IntegralGenerator
from ffcx.codegeneration.utils import dtype_to_c_type, dtype_to_scalar_dtype
from ffcx.ir.representation import IntegralIR

logger = logging.getLogger("ffcx")


def generator(ir: IntegralIR, domain: basix.CellType, options):
    """Generate C code for an integral."""
    logger.info("Generating code for integral:")
    logger.info(f"--- type: {ir.expression.integral_type}")
    logger.info(f"--- name: {ir.expression.name}")

    factory_name = f"{ir.expression.name}_{domain.name}"

    # Format declaration
    declaration = ufcx_integrals.declaration.format(factory_name=factory_name)

    # Create FFCx C backend
    backend = FFCXBackend(ir, options)

    # Configure kernel generator
    ig = IntegralGenerator(ir, backend)

    # Generate code ast for the tabulate_tensor body
    parts = ig.generate(domain)

    # Format code as string
    CF = CFormatter(options["scalar_type"])
    body = CF.c_format(parts)

    # Generate generic FFCx code snippets and add specific parts
    code = {}

    if len(ir.enabled_coefficients) > 0:
        values = ", ".join("1" if i else "0" for i in ir.enabled_coefficients)
        sizes = len(ir.enabled_coefficients)
        code["enabled_coefficients_init"] = (
            f"bool enabled_coefficients_{ir.expression.name}_{domain.name}[{sizes}] = {{{values}}};"
        )
        code["enabled_coefficients"] = f"enabled_coefficients_{ir.expression.name}_{domain.name}"
    else:
        code["enabled_coefficients_init"] = ""
        code["enabled_coefficients"] = "NULL"

    code["tabulate_tensor"] = body

    code["tabulate_tensor_float32"] = ".tabulate_tensor_float32 = NULL,"
    code["tabulate_tensor_float64"] = ".tabulate_tensor_float64 = NULL,"
    if sys.platform.startswith("win32"):
        code["tabulate_tensor_complex64"] = ""
        code["tabulate_tensor_complex128"] = ""
    else:
        code["tabulate_tensor_complex64"] = ".tabulate_tensor_complex64 = NULL,"
        code["tabulate_tensor_complex128"] = ".tabulate_tensor_complex128 = NULL,"
    np_scalar_type = np.dtype(options["scalar_type"]).name
    code[f"tabulate_tensor_{np_scalar_type}"] = (
        f".tabulate_tensor_{np_scalar_type} = tabulate_tensor_{factory_name},"
    )

    assert ir.expression.coordinate_element_hash is not None
    implementation = ufcx_integrals.factory.format(
        factory_name=factory_name,
        enabled_coefficients=code["enabled_coefficients"],
        enabled_coefficients_init=code["enabled_coefficients_init"],
        tabulate_tensor=code["tabulate_tensor"],
        needs_facet_permutations="true" if ir.expression.needs_facet_permutations else "false",
        scalar_type=dtype_to_c_type(options["scalar_type"]),
        geom_type=dtype_to_c_type(dtype_to_scalar_dtype(options["scalar_type"])),
        coordinate_element_hash=f"UINT64_C({ir.expression.coordinate_element_hash})",
        tabulate_tensor_float32=code["tabulate_tensor_float32"],
        tabulate_tensor_float64=code["tabulate_tensor_float64"],
        tabulate_tensor_complex64=code["tabulate_tensor_complex64"],
        tabulate_tensor_complex128=code["tabulate_tensor_complex128"],
        domain=int(domain),
    )

    return declaration, implementation
