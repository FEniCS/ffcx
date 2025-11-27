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
from ffcx.codegeneration.numba.numba_implementation import NumbaFormatter
from ffcx.codegeneration.utils import dtype_to_scalar_dtype
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

    d = {}
    d["name_from_uflfile"] = ir.name_from_uflfile
    d["factory_name"] = factory_name

    parts = eg.generate()

    tensor_size = 1
    for dim in ir.expression.shape:
        tensor_size *= dim

    tensor_size *= 3  # TODO: number of evaluation points - where to get?

    n_coeff = sum(coeff.ufl_element().dim for coeff in ir.expression.coefficient_offsets.keys())
    n_const = sum(
        np.prod(constant.ufl_shape, dtype=int)
        for constant in ir.expression.original_constant_offsets.keys()
    )
    n_coord_dofs = ir.expression.number_coordinate_dofs * 3
    n_entity_local_index = 2  # TODO: this is just an upper bound, harmful?
    n_quad_perm = 2 if ir.expression.needs_facet_permutations else 0

    header = f"""
    A = numba.carray(_A, ({tensor_size}))
    w = numba.carray(_w, ({n_coeff}))
    c = numba.carray(_c, ({n_const}))
    coordinate_dofs = numba.carray(_coordinate_dofs, ({n_coord_dofs}))
    entity_local_index = numba.carray(_entity_local_index, ({n_entity_local_index}))
    quadrature_permutation = numba.carray(_quadrature_permutation, ({n_quad_perm}))
    """
    F = NumbaFormatter(options["scalar_type"])
    body = F.format(parts)
    body = ["    " + line for line in body.split("\n")]
    body = "\n".join(body)

    d["tabulate_expression"] = header + body

    originals = ", ".join(str(i) for i in ir.original_coefficient_positions)
    d["original_coefficient_positions"] = f"[{originals}]"
    # n = points.size
    d["points"] = f"[{', '.join(str(p) for p in points.flatten())}]"

    shape = ", ".join(str(i) for i in ir.expression.shape)
    d["value_shape"] = f"[{shape}]"
    d["num_components"] = len(ir.expression.shape)
    d["num_coefficients"] = len(ir.expression.coefficient_numbering)
    d["num_constants"] = len(ir.constant_names)
    d["num_points"] = points.shape[0]
    d["topological_dimension"] = points.shape[1]
    d["scalar_type"] = options["scalar_type"]
    d["geom_type"] = dtype_to_scalar_dtype(options["scalar_type"])

    d["rank"] = len(ir.expression.tensor_shape)
    d["np_scalar_type"] = np.dtype(options["scalar_type"]).names

    names = ", ".join(f'"{name}"' for name in ir.coefficient_names)
    d["coefficient_names"] = f"[{names}]"
    names = ", ".join(f'"{name}"' for name in ir.constant_names)
    d["constant_names"] = f"[{names}]"

    code = []

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

    d["function_spaces_alloc"] = "\n".join(code)
    d["function_spaces"] = ""

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

    # d["function_spaces"] = "0"

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
