# Copyright (C) 2015-2021 Martin Sandve Alnæs, Michal Habera, Igor Baratta, Chris Richardson and
# Paul T. Kühner
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Generate UFCx code for an integral."""

import logging

import basix
from numpy import typing as npt

from ffcx.codegeneration.backend import FFCXBackend
from ffcx.codegeneration.common import template_keys, tensor_sizes
from ffcx.codegeneration.integral_generator import IntegralGenerator
from ffcx.codegeneration.numba import integral_template as ufcx_integrals
from ffcx.codegeneration.numba.formatter import Formatter
from ffcx.ir.representation import IntegralIR

logger = logging.getLogger("ffcx")


def generator(
    ir: IntegralIR, domain: basix.CellType, options: dict[str, int | float | npt.DTypeLike]
) -> tuple[str]:
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
    sizes = tensor_sizes(ir)
    header = f"""
    A = numba.carray(_A, ({sizes.A}))
    w = numba.carray(_w, ({sizes.w}))
    c = numba.carray(_c, ({sizes.c}))
    coordinate_dofs = numba.carray(_coordinate_dofs, ({sizes.coords}))
    entity_local_index = numba.carray(_entity_local_index, ({sizes.local_index}))
    quadrature_permutation = numba.carray(_quadrature_permutation, ({sizes.permutation}))
    """
    d["tabulate_tensor"] = header + body
    d["needs_facet_permutations"] = "True" if ir.expression.needs_facet_permutations else "False"
    d["coordinate_element_hash"] = ir.expression.coordinate_element_hash
    d["domain"] = str(int(domain))

    assert ir.expression.coordinate_element_hash is not None

    assert set(d.keys()) == template_keys(ufcx_integrals.factory)
    return (ufcx_integrals.factory.format_map(d),)
