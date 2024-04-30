# Copyright (C) 2019 Michal Habera
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import logging

from ffcx.codegeneration.numba import expressions_template
from ffcx.codegeneration.backend import FFCXBackend
from ffcx.codegeneration.expression_generator import ExpressionGenerator
from ffcx.codegeneration.numba.numba_implementation import NumbaFormatter

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

    tensor_size = 1
    for dim in ir.tensor_shape:
        tensor_size *= dim
    n_coeff = 1000
    n_const = 1000
    header = f"""
    A = numba.carray(_A, ({tensor_size}))
    w = numba.carray(_w, ({n_coeff}))
    c = numba.carray(_c, ({n_const}))
    coordinate_dofs = numba.carray(_coordinate_dofs, (1000))
    entity_local_index = numba.carray(_entity_local_index, (1000))
    quadrature_permutation = numba.carray(_quadrature_permutation, (1000))
    """
    F = NumbaFormatter(options["scalar_type"])
    body = F.c_format(parts)
    body = ["    " + line for line in body.split("\n")]
    body = "\n".join(body)

    d["tabulate_expression"] = header + body

    originals = ", ".join(str(i) for i in ir.original_coefficient_positions)
    d["original_coefficient_positions"] = f"[{originals}]"
    points = ", ".join(str(p) for p in ir.points.flatten())
    n = ir.points.size
    d["points"] = f"[{points}]"

    shape = ", ".join(str(i) for i in ir.expression_shape)
    d["value_shape"] = f"[{shape}]"
    d["num_components"] = len(ir.expression_shape)
    d["num_coefficients"] = len(ir.coefficient_numbering)
    d["num_constants"] = len(ir.constant_names)
    d["num_points"] = ir.points.shape[0]
    d["topological_dimension"] = ir.points.shape[1]
    d["scalar_type"] = options["scalar_type"]
    
    d["rank"] = len(ir.tensor_shape)

    names = ", ".join(f'"{name}"' for name in ir.coefficient_names)
    d["coefficient_names"] = f"[{names}]"
    names = ", ".join(f'"{name}"' for name in ir.constant_names)
    d["constant_names"] = f"[{names}]"

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

    d["function_spaces"] = "0"

    # Check that no keys are redundant or have been missed
    # from string import Formatter

    # fields = [
    #     fname
    #     for _, fname, _, _ in Formatter().parse(expressions_template.factory)
    #     if fname
    # ]
    # assert set(fields) == set(
    #     d.keys()
    # ), "Mismatch between keys in template and in formatting dict"

    # Format implementation code
    implementation = expressions_template.factory.format_map(d)

    return declaration, implementation
