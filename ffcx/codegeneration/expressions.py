# Copyright (C) 2019 Michal Habera
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

from ffcx.codegeneration import expressions_template
from ffcx.codegeneration.expressions_generator import generate_expression_code


def generator(ir, parameters):
    """Generate UFC code for an expression."""
    factory_name = ir.name

    # Format declaration
    declaration = expressions_template.declaration.format(factory_name=factory_name)

    # Generate code
    code = generate_expression_code(ir, parameters)

    # Format implementation code
    implementation = expressions_template.factory.format(
        factory_name=factory_name,
        tabulate_expression=code["tabulate_expression"],
        original_coefficient_positions=code["original_coefficient_positions"],
        num_coefficients=len(ir.coefficient_numbering),
        num_points=ir.points.shape[0],
        topological_dimension=ir.points.shape[1],
        num_components=len(ir.expression_shape),
        points=code["points"],
        value_shape=code["value_shape"])

    return declaration, implementation
