# -*- coding: utf-8 -*-
# Copyright (C) 2015-2017 Martin Sandve Alnæs
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# You should have received a copy of the GNU Lesser General Public License
# along with UFLACS. If not, see <http://www.gnu.org/licenses/>.

from ffc.codegeneration import integrals_template as ufc_integrals


def generator(ir, parameters):
    """Generate UFC code for an integral"""
    factory_name = ir.classname
    integral_type = ir.integral_type

    # Format declaration
    if integral_type == "custom":
        declaration = ufc_integrals.custom_declaration.format(factory_name=factory_name)
    else:
        declaration = ufc_integrals.declaration.format(factory_name=factory_name)

    if ir.representation == "uflacs":
        from ffc.codegeneration.uflacsgenerator import generate_integral_code
    elif ir.representation == "tsfc":
        from ffc.codegeneration.tsfcgenerator import generate_integral_code
    else:
        raise RuntimeError("Unknown representation: {}".format(ir.representation))

    # Generate code
    code = generate_integral_code(ir, parameters)

    # Format tabulate tensor body
    tabulate_tensor_declaration = ufc_integrals.tabulate_implementation[
        integral_type]
    tabulate_tensor_fn = tabulate_tensor_declaration.format(
        factory_name=factory_name, tabulate_tensor=code["tabulate_tensor"])

    # Format implementation code

    if integral_type == "custom":
        implementation = ufc_integrals.custom_factory.format(
            factory_name=factory_name,
            enabled_coefficients=code["enabled_coefficients"],
            tabulate_tensor=tabulate_tensor_fn)
    else:
        implementation = ufc_integrals.factory.format(
            factory_name=factory_name,
            enabled_coefficients=code["enabled_coefficients"],
            tabulate_tensor=tabulate_tensor_fn)

    return declaration, implementation
