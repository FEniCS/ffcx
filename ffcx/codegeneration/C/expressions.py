# Copyright (C) 2019 Michal Habera
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Generate UFC code for an expression."""

from __future__ import annotations

import logging

import numpy as np

from ffcx.codegeneration.backend import FFCXBackend
from ffcx.codegeneration.C import expressions_template
from ffcx.codegeneration.C.c_implementation import CFormatter
from ffcx.codegeneration.expression_generator import ExpressionGenerator
from ffcx.codegeneration.utils import dtype_to_c_type, dtype_to_scalar_dtype
from ffcx.ir.representation import ExpressionIR

logger = logging.getLogger("ffcx")


def generator(ir: ExpressionIR, options):
    """Generate UFC code for an expression."""
    logger.info("Generating code for expression:")
    assert len(ir.expression.integrand) == 1, "Expressions only support single quadrature rule"
    points = next(iter(ir.expression.integrand)).points
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

    CF = CFormatter(options["scalar_type"])
    d["tabulate_expression"] = CF.c_format(parts)

    if len(ir.original_coefficient_positions) > 0:
        d["original_coefficient_positions"] = f"original_coefficient_positions_{factory_name}"
        values = ", ".join(str(i) for i in ir.original_coefficient_positions)
        sizes = len(ir.original_coefficient_positions)
        d["original_coefficient_positions_init"] = (
            f"static int original_coefficient_positions_{factory_name}[{sizes}] = {{{values}}};"
        )
    else:
        d["original_coefficient_positions"] = "NULL"
        d["original_coefficient_positions_init"] = ""

    values = ", ".join(str(p) for p in points.flatten())
    sizes = points.size
    d["points_init"] = f"static double points_{factory_name}[{sizes}] = {{{values}}};"
    d["points"] = f"points_{factory_name}"

    if len(ir.expression.shape) > 0:
        values = ", ".join(str(i) for i in ir.expression.shape)
        sizes = len(ir.expression.shape)
        d["value_shape_init"] = f"static int value_shape_{factory_name}[{sizes}] = {{{values}}};"
        d["value_shape"] = f"value_shape_{factory_name}"
    else:
        d["value_shape_init"] = ""
        d["value_shape"] = "NULL"
    d["num_components"] = len(ir.expression.shape)
    d["num_coefficients"] = len(ir.expression.coefficient_numbering)
    d["num_constants"] = len(ir.constant_names)
    d["num_points"] = points.shape[0]
    d["entity_dimension"] = points.shape[1]
    d["scalar_type"] = dtype_to_c_type(options["scalar_type"])
    d["geom_type"] = dtype_to_c_type(dtype_to_scalar_dtype(options["scalar_type"]))
    d["np_scalar_type"] = np.dtype(options["scalar_type"]).name

    d["rank"] = len(ir.expression.tensor_shape)

    if len(ir.coefficient_names) > 0:
        values = ", ".join(f'"{name}"' for name in ir.coefficient_names)
        sizes = len(ir.coefficient_names)
        d["coefficient_names_init"] = (
            f"static const char* coefficient_names_{factory_name}[{sizes}] = {{{values}}};"
        )
        d["coefficient_names"] = f"coefficient_names_{factory_name}"
    else:
        d["coefficient_names_init"] = ""
        d["coefficient_names"] = "NULL"

    if len(ir.constant_names) > 0:
        values = ", ".join(f'"{name}"' for name in ir.constant_names)
        sizes = len(ir.constant_names)
        d["constant_names_init"] = (
            f"static const char* constant_names_{factory_name}[{sizes}] = {{{values}}};"
        )
        d["constant_names"] = f"constant_names_{factory_name}"
    else:
        d["constant_names_init"] = ""
        d["constant_names"] = "NULL"

    # Check that no keys are redundant or have been missed
    from string import Formatter

    fields = [fname for _, fname, _, _ in Formatter().parse(expressions_template.factory) if fname]
    assert set(fields) == set(d.keys()), "Mismatch between keys in template and in formatting dict"

    # Format implementation code
    implementation = expressions_template.factory.format_map(d)

    return declaration, implementation
