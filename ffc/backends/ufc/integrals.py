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


class Snippets:
    "Each function maps to a keyword in the template.."

    def enabled_coefficients(self, L, ir):
        # FIXME: Needs updating for interface change to plain C.
        enabled_coefficients = ir["enabled_coefficients"]
        initializer_list = ", ".join(
            "true" if enabled else "false" for enabled in enabled_coefficients)
        code = L.StatementList([
            # Cheating a bit with verbatim:
            L.VerbatimStatement("static const std::vector<bool> enabled({%s});"
                                % initializer_list),
            L.Return(L.Symbol("enabled")),
        ])
        return code

    def tabulate_tensor(self, L, ir):
        # FIXME: This is where the current ffc backend code generation should be injected,
        #        however that will require using pick_representation in here?
        tt = ir["tabulate_tensor"]
        code = "code generated from %s" % tt
        return code

    def tabulate_tensor_comment(self, L, ir):
        "Generate comment for tabulate_tensor."

        r = ir["representation"]
        integrals_metadata = ir["integrals_metadata"]
        integral_metadata = ir["integral_metadata"]

        lines = [
            "This function was generated using '%s' representation" % r,
            "with the following integrals metadata:", ""
        ]
        lines += [("  " + l) for l in dstr(integrals_metadata).split("\n")]
        lines += [""]
        for i, metadata in enumerate(integral_metadata):
            lines += [
                "",
                "and the following integral %d metadata:" % i,
                "",
            ]
            lines += [("  " + l) for l in dstr(metadata).split("\n")][:-1]
            lines += [""]

        return [L.Comment(line) for line in lines]


class ufc_cell_integral(Snippets):
    def __init__(self):
        Snippets.__init__(self, "cell")


class ufc_exterior_facet_integral(Snippets):
    def __init__(self):
        Snippets.__init__(self, "exterior_facet")


class ufc_interior_facet_integral(Snippets):
    def __init__(self):
        Snippets.__init__(self, "interior_facet")


class ufc_vertex_integral(Snippets):
    def __init__(self):
        Snippets.__init__(self, "vertex")


class ufc_custom_integral(Snippets):
    def __init__(self):
        Snippets.__init__(self, "custom")

    def num_cells(self, L, ir):
        value = ir["num_cells"]
        return L.Return(L.LiteralInt(value))


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
