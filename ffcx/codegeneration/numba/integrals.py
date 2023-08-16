# Copyright (C) 2015-2021 Martin Sandve Alnæs, Michal Habera, Igor Baratta
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import logging

from ffcx.codegeneration.integral_generator import IntegralGenerator
from ffcx.codegeneration.numba import integrals_template as ufcx_integrals
from ffcx.codegeneration.backend import FFCXBackend
from ffcx.codegeneration.utils import cdtype_to_numpy, scalar_to_value_type
from ffcx.codegeneration.numba.numba_implementation import NumbaFormatter


logger = logging.getLogger("ffcx")


def generator(ir, options):
    logger.info("Generating code for integral:")
    logger.info(f"--- type: {ir.integral_type}")
    logger.info(f"--- name: {ir.name}")

    """Generate code for an integral."""
    factory_name = ir.name

    # Create FFCx backend
    backend = FFCXBackend(ir, options)

    # Configure kernel generator
    ig = IntegralGenerator(ir, backend)

    # Generate code AST for the tabulate_tensor body
    parts = ig.generate()

    # Format code as string
    F = NumbaFormatter(options["scalar_type"])
    body = F.c_format(parts)
    body = ["    " + line for line in body.split("\n")]
    body = "\n".join(body)

    # Generate generic FFCx code snippets and add specific parts
    code = {}
    code["class_type"] = ir.integral_type + "_integral"
    code["name"] = ir.name

    vals = ", ".join("1" if i else "0" for i in ir.enabled_coefficients)
    code["enabled_coefficients"] = f"[{vals}]"

    tensor_size = ir.tensor_shape[0]
    if len(ir.tensor_shape) == 2:
        tensor_size *= ir.tensor_shape[1]
    n_coeff = 1000
    n_const = 1000
    header = f"    A = numba.carray(_A, ({tensor_size}))\n"
    header += f"    w = numba.carray(_w, ({n_coeff}))\n"
    header += f"    c = numba.carray(_c, ({n_const}))\n"
    header += "    coordinate_dofs = numba.carray(_coordinate_dofs, (1000))\n"
    code["tabulate_tensor"] = header + body

    implementation = ufcx_integrals.factory.format(
        factory_name=factory_name,
        enabled_coefficients=code["enabled_coefficients"],
        tabulate_tensor=code["tabulate_tensor"],
        needs_facet_permutations="True" if ir.needs_facet_permutations else "False",
        scalar_type=options["scalar_type"],
        geom_type=scalar_to_value_type(options["scalar_type"]),
        np_scalar_type=cdtype_to_numpy(options["scalar_type"]),
        coordinate_element=ir.coordinate_element,
    )

    return "", implementation
