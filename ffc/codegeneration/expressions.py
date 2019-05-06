# -*- coding: utf-8 -*-
# Copyright (C) 2019 Michal Habera
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

from ffc.codegeneration import expressions_template
from ffc.codegeneration.uflacsgenerator import generate_expression_code

def generator(ir, parameters):
    """Generate UFC code for an expression"""
    factory_name = ir.classname

    # Format declaration
    declaration = expressions_template.declaration.format(factory_name=factory_name)

    # Generate code
    code = generate_expression_code(ir, parameters)

    # Format tabulate tensor body
    tabulate_expression_declaration = expressions_template.tabulate_implementation
    tabulate_expression_fn = tabulate_expression_declaration.format(
        factory_name=factory_name, tabulate_expression=code["tabulate_expression"])

    # Format implementation code
    implementation = expressions_template.factory.format(
        factory_name=factory_name,
        tabulate_expression=tabulate_expression_fn)

    return declaration, implementation
