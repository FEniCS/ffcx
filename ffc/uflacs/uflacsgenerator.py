# -*- coding: utf-8 -*-
# Copyright (C) 2013-2016 Martin Sandve Alnæs
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.

"""Controlling algorithm for building the tabulate_tensor
source structure from factorized representation."""

from ffc.log import info
from ffc.representationutils import initialize_integral_code

from ffc.uflacs.backends.ffc.backend import FFCBackend
from ffc.uflacs.generation.integralgenerator import IntegralGenerator
from ffc.uflacs.language.format_lines import format_indented_lines


def generate_integral_code(ir, prefix, parameters):
    "Generate code for integral from intermediate representation."

    info("Generating code from ffc.uflacs representation")

    # Generate generic ffc code snippets
    code = initialize_integral_code(ir, prefix, parameters)

    # Generate tabulate_tensor body using uflacs algorithms
    uflacs_code = generate_tabulate_tensor_code(ir, prefix, parameters)

    code["tabulate_tensor"] = uflacs_code["tabulate_tensor"]

    # TODO: Use code generation utils here for consistency
    if ir.get("num_cells") is not None:
        code["num_cells"] = "  return %d;" % (ir["num_cells"],)

    code["additional_includes_set"] = set()
    code["additional_includes_set"].update(ir.get("additional_includes_set",()))
    code["additional_includes_set"].update(uflacs_code["additional_includes_set"])

    return code


def generate_tabulate_tensor_code(ir, prefix, parameters):

    # Create FFC C++ backend
    backend = FFCBackend(ir, parameters)

    # Create code generator for integral body
    ig = IntegralGenerator(ir, backend)

    # Generate code ast for the tabulate_tensor body
    parts = ig.generate()

    # Format code AST as one string
    body = format_indented_lines(parts.cs_format(), 1)

    # Fetch includes
    includes = set(ig.get_includes())

    # Format uflacs specific code structures into a single
    # string and place in dict before returning to ffc
    code = {
        "tabulate_tensor": body,
        "additional_includes_set": includes,
    }

    return code
