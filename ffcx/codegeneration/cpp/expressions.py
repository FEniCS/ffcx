# Copyright (C) 2019 Michal Habera
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Generate code for an expression."""

import logging
import string

import numpy as np

from ffcx.codegeneration.backend import FFCXBackend
from ffcx.codegeneration.C import expressions_template
from ffcx.codegeneration.C.implementation import Formatter
from ffcx.codegeneration.expression_generator import ExpressionGenerator
from ffcx.codegeneration.utils import dtype_to_c_type, dtype_to_scalar_dtype

logger = logging.getLogger("ffcx")


def generator(ir, options):
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

    d = {}
    d["name_from_uflfile"] = ir.name_from_uflfile
    d["factory_name"] = factory_name

    parts = eg.generate()

    CF = Formatter(options["scalar_type"])
    d["tabulate_expression"] = CF.format(parts)

    if len(ir.original_coefficient_positions) > 0:
        d["original_coefficient_positions"] = f"original_coefficient_positions_{factory_name}"
        sizes = len(ir.original_coefficient_positions)
        values = ", ".join(str(i) for i in ir.original_coefficient_positions)
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

    # TODO: make cpp
    d["coordinate_element_hash"] = f"UINT64_C({ir.expression.coordinate_element_hash})"

    # FIXME: Should be handled differently, revise how
    # ufcx_function_space is generated (also for ufcx_form)
    # for name, (element, dofmap, cmap_family, cmap_degree) in ir.function_spaces.items():
    #     code += [f"static ufcx_function_space function_space_{name}_{ir.name_from_uflfile} ="]
    #     code += ["{"]
    #     code += [f".finite_element = &{element},"]
    #     code += [f".dofmap = &{dofmap},"]
    #     code += [f'.geometry_family = "{cmap_family}",']
    #     code += [f".geometry_degree = {cmap_degree}"]
    #     code += ["};"]

    # d["function_spaces_alloc"] = "\n".join(code)
    # d["function_spaces"] = ""

    # if len(ir.function_spaces) > 0:
    #     d["function_spaces"] = f"function_spaces_{ir.name}"
    #     fs_list = ", ".join(
    #         f"&function_space_{name}_{ir.name_from_uflfile}"
    #         for (name, _) in ir.function_spaces.items()
    #     )
    #     n = len(ir.function_spaces.items())
    #     d["function_spaces_init"] = (
    #         f"ufcx_function_space* function_spaces_{ir.name}[{n}] = {{{fs_list}}};"
    #     )
    # else:
    #     d["function_spaces"] = "NULL"
    #     d["function_spaces_init"] = ""

    # Check that no keys are redundant or have been missed
    fields = [
        fname for _, fname, _, _ in string.Formatter().parse(expressions_template.factory) if fname
    ]
    assert set(fields) == set(d.keys()), "Mismatch between keys in template and in formatting dict"

    # Format implementation code
    implementation = expressions_template.factory.format_map(d)

    return declaration, implementation
