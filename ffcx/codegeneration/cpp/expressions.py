# Copyright (C) 2019 Michal Habera
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Generate code for an expression."""

import logging

from ffcx.codegeneration.backend import FFCXBackend
from ffcx.codegeneration.C import expressions_template
from ffcx.codegeneration.C.c_implementation import CFormatter
from ffcx.codegeneration.expression_generator import ExpressionGenerator

logger = logging.getLogger("ffcx")


def generator(ir, options):
    """Generate UFC code for an expression."""
    logger.info("Generating code for expression:")
    logger.info(f"--- points: {ir.points}")
    logger.info(f"--- name: {ir.name}")

    factory_name = ir.name

    # Format declaration
    declaration = expressions_template.declaration.format(
        factory_name=factory_name, name_from_uflfile=ir.name_from_uflfile
    )

    backend = FFCXBackend(ir, options)
    eg = ExpressionGenerator(ir, backend)

    d = {}
    d["name_from_uflfile"] = ir.name_from_uflfile
    d["factory_name"] = ir.name

    parts = eg.generate()

    CF = CFormatter(options["scalar_type"])
    d["tabulate_expression"] = CF.c_format(parts)

    if len(ir.original_coefficient_positions) > 0:
        d[
            "original_coefficient_positions"
        ] = f"original_coefficient_positions_{ir.name}"
        n = len(ir.original_coefficient_positions)
        originals = ", ".join(str(i) for i in ir.original_coefficient_positions)
        d[
            "original_coefficient_positions_init"
        ] = f"static int original_coefficient_positions_{ir.name}[{n}] = {{{originals}}};"

    else:
        d["original_coefficient_positions"] = "NULL"
        d["original_coefficient_positions_init"] = ""

    points = ", ".join(str(p) for p in ir.points.flatten())
    n = ir.points.size
    d["points_init"] = f"static double points_{ir.name}[{n}] = {{{points}}};"
    d["points"] = f"points_{ir.name}"

    if len(ir.expression_shape) > 0:
        n = len(ir.expression_shape)
        shape = ", ".join(str(i) for i in ir.expression_shape)
        d["value_shape_init"] = f"static int value_shape_{ir.name}[{n}] = {{{shape}}};"
        d["value_shape"] = f"value_shape_{ir.name}"
    else:
        d["value_shape_init"] = ""
        d["value_shape"] = "NULL"

    d["num_components"] = len(ir.expression_shape)
    d["num_coefficients"] = len(ir.coefficient_numbering)
    d["num_constants"] = len(ir.constant_names)
    d["num_points"] = ir.points.shape[0]
    d["topological_dimension"] = ir.points.shape[1]
    d["rank"] = len(ir.tensor_shape)

    if len(ir.coefficient_names) > 0:
        names = ", ".join(f'"{name}"' for name in ir.coefficient_names)
        n = len(ir.coefficient_names)
        d[
            "coefficient_names_init"
        ] = f"static const char* coefficient_names_{ir.name}[{n}] = {{{names}}};"

        d["coefficient_names"] = f"coefficient_names_{ir.name}"
    else:
        d["coefficient_names_init"] = ""
        d["coefficient_names"] = "NULL"

    if len(ir.constant_names) > 0:
        names = ", ".join(f'"{name}"' for name in ir.constant_names)
        n = len(ir.constant_names)
        d[
            "constant_names_init"
        ] = f"static const char* constant_names_{ir.name}[{n}] = {{{names}}};"
        d["constant_names"] = f"constant_names_{ir.name}"
    else:
        d["constant_names_init"] = ""
        d["constant_names"] = "NULL"

    code = []

    # FIXME: Should be handled differently, revise how
    # ufcx_function_space is generated (also for ufcx_form)
    for name, (element, dofmap, cmap_family, cmap_degree) in ir.function_spaces.items():
        code += [
            f"static ufcx_function_space function_space_{name}_{ir.name_from_uflfile} ="
        ]
        code += ["{"]
        code += [f".finite_element = &{element},"]
        code += [f".dofmap = &{dofmap},"]
        code += [f'.geometry_family = "{cmap_family}",']
        code += [f".geometry_degree = {cmap_degree}"]
        code += ["};"]

    d["function_spaces_alloc"] = "\n".join(code)
    d["function_spaces"] = ""

    if len(ir.function_spaces) > 0:
        d["function_spaces"] = f"function_spaces_{ir.name}"
        fs_list = ", ".join(
            f"&function_space_{name}_{ir.name_from_uflfile}"
            for (name, _) in ir.function_spaces.items()
        )
        n = len(ir.function_spaces.items())
        d[
            "function_spaces_init"
        ] = f"ufcx_function_space* function_spaces_{ir.name}[{n}] = {{{fs_list}}};"
    else:
        d["function_spaces"] = "NULL"
        d["function_spaces_init"] = ""

    # Check that no keys are redundant or have been missed
    from string import Formatter

    fields = [
        fname
        for _, fname, _, _ in Formatter().parse(expressions_template.factory)
        if fname
    ]
    assert set(fields) == set(
        d.keys()
    ), "Mismatch between keys in template and in formatting dict"

    # Format implementation code
    implementation = expressions_template.factory.format_map(d)

    return declaration, implementation
