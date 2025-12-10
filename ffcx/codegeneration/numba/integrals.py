# Copyright (C) 2015-2021 Martin Sandve Alnæs, Michal Habera, Igor Baratta, Chris Richardson and
# Paul T. Kühner
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Generate UFC code for an integral."""

import logging
import string

import basix
import numpy as np
from numpy import typing as npt

from ffcx.codegeneration.backend import FFCXBackend
from ffcx.codegeneration.integral_generator import IntegralGenerator
from ffcx.codegeneration.numba import integrals_template as ufcx_integrals
from ffcx.codegeneration.numba.implementation import Formatter
from ffcx.ir.representation import IntegralIR

logger = logging.getLogger("ffcx")


def generator(
    ir: IntegralIR, domain: basix.CellType, options: dict[str, int | float | npt.DTypeLike]
) -> tuple[str, str]:
    """Generate numba code for an integral.

    Args:
        ir: IR of the integral
        domain: basix cell type
        options: dict of kernel generation options

    Returns:
        Tuple of declaration (header) and implementation (source) strings.

    Note:
        numba backend only provides a declaration. Implementation string will always be empty.

    """
    logger.info("Generating code for integral:")
    logger.info(f"--- type: {ir.expression.integral_type}")
    logger.info(f"--- name: {ir.expression.name}")

    factory_name = f"{ir.expression.name}_{domain.name}"

    # Format declaration
    declaration = ""

    # Create FFCx backend
    backend = FFCXBackend(ir, options)

    # Configure kernel generator
    ig = IntegralGenerator(ir, backend)

    # Generate code AST for the tabulate_tensor body
    parts = ig.generate(domain)

    # Format code as string
    format = Formatter(options["scalar_type"])  # type: ignore
    body = format(parts)
    body = "\n".join(["    " + line for line in body.split("\n")])

    # Generate generic FFCx code snippets and add specific parts
    d: dict[str, str] = {}

    d["factory_name"] = factory_name

    # TODO: enabled_coefficients_init - required?
    vals = ", ".join("1" if i else "0" for i in ir.enabled_coefficients)
    d["enabled_coefficients"] = f"[{vals}]"

    # tabulate_tensor
    # Note: In contrast to the C implementation we actually need to provide/compute the sizes of the
    #       array.
    size_A = np.prod(ir.expression.tensor_shape)
    size_w = sum(coeff.ufl_element().dim for coeff in ir.expression.coefficient_offsets.keys())
    size_c = sum(
        np.prod(constant.ufl_shape, dtype=int)
        for constant in ir.expression.original_constant_offsets.keys()
    )
    size_coords = ir.expression.number_coordinate_dofs * 3
    size_local_index = 2  # TODO: this is just an upper bound
    size_permutation = 2 if ir.expression.needs_facet_permutations else 0

    header = f"""
    A = numba.carray(_A, ({size_A}))
    w = numba.carray(_w, ({size_w}))
    c = numba.carray(_c, ({size_c}))
    coordinate_dofs = numba.carray(_coordinate_dofs, ({size_coords}))
    entity_local_index = numba.carray(_entity_local_index, ({size_local_index}))
    quadrature_permutation = numba.carray(_quadrature_permutation, ({size_permutation}))
    """
    d["tabulate_tensor"] = header + body
    d["needs_facet_permutations"] = "True" if ir.expression.needs_facet_permutations else "False"
    d["coordinate_element_hash"] = ir.expression.coordinate_element_hash
    d["domain"] = str(int(domain))

    assert ir.expression.coordinate_element_hash is not None
    implementation = ufcx_integrals.factory.format_map(d)

    # Check that no keys are redundant or have been missed
    fields = [fname for _, fname, _, _ in string.Formatter().parse(ufcx_integrals.factory) if fname]
    assert set(fields) == set(d.keys()), "Mismatch between keys in template and in formatting dict"

    return declaration, implementation
