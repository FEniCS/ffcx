# Copyright (C) 2019-2025 Michal Habera, Chris Richardson and Paul T. Kühner
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Generate UFC code for an expression."""

import logging

import numpy as np
import numpy.typing as npt

from ffcx.codegeneration.backend import FFCXBackend
from ffcx.codegeneration.common import template_keys, tensor_sizes
from ffcx.codegeneration.expression_generator import ExpressionGenerator
from ffcx.codegeneration.numba import expressions_template
from ffcx.codegeneration.numba.formatter import Formatter
from ffcx.ir.representation import ExpressionIR

logger = logging.getLogger("ffcx")


def generator(ir: ExpressionIR, options: dict[str, int | float | npt.DTypeLike]) -> tuple[str, str]:
    """Generate UFC code for an expression."""
    logger.info("Generating code for expression:")
    assert len(ir.expression.integrand) == 1, "Expressions only support single quadrature rule"
    points = next(iter(ir.expression.integrand))[1].points
    logger.info(f"--- points: {points}")
    factory_name = ir.expression.name
    logger.info(f"--- name: {factory_name}")

    # Format declaration
    declaration = expressions_template.declaration.format(
        factory_name=factory_name, name_from_uflfile=ir.name_from_uflfile
    )

    backend = FFCXBackend(ir, options)
    eg = ExpressionGenerator(ir, backend)

    d: dict[str, str | int] = {}
    d["name_from_uflfile"] = ir.name_from_uflfile
    d["factory_name"] = factory_name
    parts = eg.generate()

    # tabulate_expression
    sizes = tensor_sizes(ir)
    header = f"""
    A = numba.carray(_A, ({sizes.A}))
    w = numba.carray(_w, ({sizes.w}))
    c = numba.carray(_c, ({sizes.c}))
    coordinate_dofs = numba.carray(_coordinate_dofs, ({sizes.coords}))
    entity_local_index = numba.carray(_entity_local_index, ({sizes.local_index}))
    quadrature_permutation = numba.carray(_quadrature_permutation, ({sizes.permutation}))
    """
    format = Formatter(options["scalar_type"])  # type: ignore
    body = format(parts)
    body = "\n".join(["    " + line for line in body.split("\n")])

    d["tabulate_expression"] = header + body

    # TODO: original_coefficient_positions_init
    originals = ", ".join(str(i) for i in ir.original_coefficient_positions)
    d["original_coefficient_positions"] = f"[{originals}]"

    # TODO: points_init
    d["points"] = f"[{', '.join(str(p) for p in points.flatten())}]"

    # TODO: value_shape_init
    shape = ", ".join(str(i) for i in ir.expression.shape)
    d["value_shape"] = f"[{shape}]"
    d["num_components"] = int(np.prod(ir.expression.shape, dtype=np.int32))
    d["num_coefficients"] = len(ir.expression.coefficient_numbering)
    d["num_constants"] = len(ir.constant_names)
    d["num_points"] = points.shape[0]
    d["entity_dimension"] = points.shape[1]

    d["rank"] = len(ir.expression.tensor_shape)

    # TODO: coefficient_names_init
    names = ", ".join(f'"{name}"' for name in ir.coefficient_names)
    d["coefficient_names"] = f"[{names}]"

    # TODO: constant_names_init
    names = ", ".join(f'"{name}"' for name in ir.constant_names)
    d["constant_names"] = f"[{names}]"

    d["coordinate_element_hash"] = ir.expression.coordinate_element_hash

    # Format implementation code
    assert set(d.keys()) == template_keys(expressions_template.factory)
    implementation = expressions_template.factory.format_map(d)

    return declaration, implementation
