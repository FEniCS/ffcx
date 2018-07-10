# -*- coding: utf-8 -*-
# Copyright (C) 2015-2017 Martin Sandve Alnæs
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# You should have received a copy of the GNU Lesser General Public License
# along with UFLACS. If not, see <http://www.gnu.org/licenses/>.

from ufl.utils.formatting import dstr

from ffc.representation import pick_representation
from ffc.backends.ufc import integrals_template as ufc_integrals


def ufc_integral_generator(ir, parameters):
    """Generate UFC code for an integral"""
    factory_name = ir["classname"]
    integral_type = ir["integral_type"]

    # Format declaration
    declaration = ufc_integrals.declaration.format(
        type=integral_type, factory_name=factory_name)

    # Select representation
    r = pick_representation(ir["representation"])

    # Generate code
    # TODO: Drop prefix argument and get from ir:
    code = r.generate_integral_code(ir, ir["prefix"], parameters)

    # Hack for benchmarking overhead in assembler with empty
    # tabulate_tensor
    if parameters["generate_dummy_tabulate_tensor"]:
        code["tabulate_tensor"] = ""

    # Format tabulate tensor body
    tabulate_tensor_declaration = ufc_integrals.tabulate_implementation[
        integral_type]
    tabulate_tensor_fn = tabulate_tensor_declaration.format(
        factory_name=factory_name, tabulate_tensor=code["tabulate_tensor"])

    # Format implementation code
    implementation = ufc_integrals.factory.format(
        type=integral_type,
        factory_name=factory_name,
        enabled_coefficients=code["enabled_coefficients"],
        tabulate_tensor=tabulate_tensor_fn)

    return declaration, implementation
