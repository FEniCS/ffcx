# Copyright (C) 2015-2021 Martin Sandve Alnæs, Michal Habera, Igor Baratta
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Generate UFC code for an integral."""

import logging

from ffcx.codegeneration.backend import FFCXBackend
from ffcx.codegeneration.integral_generator import IntegralGenerator
from ffcx.codegeneration.numba import integrals_template as ufcx_integrals
from ffcx.codegeneration.numba.numba_implementation import NumbaFormatter

logger = logging.getLogger("ffcx")


def generator(ir, options):
    """Integral generator."""
    logger.info("Generating code for integral:")
    logger.info(f"--- type: {ir.expression.integral_type}")
    logger.info(f"--- name: {ir.expression.name}")

    """Generate code for an integral."""
    factory_name = ir.expression.name

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
    code["class_type"] = ir.expression.integral_type + "_integral"
    code["name"] = ir.expression.name

    vals = ", ".join("1" if i else "0" for i in ir.enabled_coefficients)
    code["enabled_coefficients"] = f"[{vals}]"

    tensor_size = 1
    for dim in ir.expression.tensor_shape:
        tensor_size *= dim
    n_coeff = len(ir.enabled_coefficients)
    n_const = len(ir.expression.original_constant_offsets)

    header = f"""
    A = numba.carray(_A, ({tensor_size}))
    w = numba.carray(_w, ({n_coeff}))
    c = numba.carray(_c, ({n_const}))
    coordinate_dofs = numba.carray(_coordinate_dofs, (1000))
    entity_local_index = numba.carray(_entity_local_index, (1000))
    quadrature_permutation = numba.carray(_quadrature_permutation, (1000))
    """
    code["tabulate_tensor"] = header + body

    implementation = ufcx_integrals.factory.format(
        factory_name=factory_name,
        enabled_coefficients=code["enabled_coefficients"],
        tabulate_tensor=code["tabulate_tensor"],
        needs_facet_permutations="True" if ir.expression.needs_facet_permutations else "False",
        scalar_type=options["scalar_type"],
        coordinate_element=ir.coordinate_element_hash,
    )

    return "", implementation
