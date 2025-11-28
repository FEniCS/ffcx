# Copyright (C) 2019-2025 Michal Habera, Chris Richardson and Paul T. Kühner
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Generate UFC code for an expression."""

import logging

import numpy as np

from ffcx.codegeneration.backend import FFCXBackend
from ffcx.codegeneration.expression_generator import ExpressionGenerator
from ffcx.codegeneration.numba import expressions_template
from ffcx.codegeneration.numba.implementation import Formatter

# from ffcx.codegeneration.utils import dtype_to_scalar_dtype
from ffcx.ir.representation import ExpressionIR

logger = logging.getLogger("ffcx")


def generator(ir: ExpressionIR, options):
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
    size_w = sum(coeff.ufl_element().dim for coeff in ir.expression.coefficient_offsets.keys())
    size_c = sum(
        np.prod(constant.ufl_shape, dtype=int)
        for constant in ir.expression.original_constant_offsets.keys()
    )
    size_coords = ir.expression.number_coordinate_dofs * 3
    size_local_index = 2  # TODO: this is just an upper bound, harmful?
    size_permutation = 2 if ir.expression.needs_facet_permutations else 0

    header = f"""
    A = numba.carray(_A, ({points.size}))
    w = numba.carray(_w, ({size_w}))
    c = numba.carray(_c, ({size_c}))
    coordinate_dofs = numba.carray(_coordinate_dofs, ({size_coords}))
    entity_local_index = numba.carray(_entity_local_index, ({size_local_index}))
    quadrature_permutation = numba.carray(_quadrature_permutation, ({size_permutation}))
    """
    F = Formatter(options["scalar_type"])
    body = F.format(parts)
    body = ["    " + line for line in body.split("\n")]
    body = "\n".join(body)

    d["tabulate_expression"] = header + body

    # TODO: original_coefficient_positions_init
    originals = ", ".join(str(i) for i in ir.original_coefficient_positions)
    d["original_coefficient_positions"] = f"[{originals}]"

    # TODO: points_init
    d["points"] = f"[{', '.join(str(p) for p in points.flatten())}]"

    # TODO: value_shape_init
    shape = ", ".join(str(i) for i in ir.expression.shape)
    d["value_shape"] = f"[{shape}]"
    d["num_components"] = len(ir.expression.shape)
    d["num_coefficients"] = len(ir.expression.coefficient_numbering)
    d["num_constants"] = len(ir.constant_names)
    d["num_points"] = points.shape[0]
    d["entity_dimension"] = points.shape[1]
    d["scalar_type"] = options["scalar_type"]
    # d["geom_type"] = dtype_to_scalar_dtype(options["scalar_type"])
    d["np_scalar_type"] = np.dtype(options["scalar_type"]).names

    d["rank"] = len(ir.expression.tensor_shape)

    # TODO: coefficient_names_init
    names = ", ".join(f'"{name}"' for name in ir.coefficient_names)
    d["coefficient_names"] = f"[{names}]"

    # TODO: constant_names_init
    names = ", ".join(f'"{name}"' for name in ir.constant_names)
    d["constant_names"] = f"[{names}]"

    # TODO: coordinate_element_hash

    # Check that no keys are redundant or have been missed

    # fields = [
    #     fname for _, fname, _, _ in string.Formatter().parse(expressions_template.factory) if
    # fname
    # ]
    # assert set(fields) == set(d.keys()), "Mismatch between keys in template and in formatting
    # dict"

    # Format implementation code
    implementation = expressions_template.factory.format_map(d)

    return declaration, implementation
