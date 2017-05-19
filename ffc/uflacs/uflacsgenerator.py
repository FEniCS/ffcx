# -*- coding: utf-8 -*-
# Copyright (C) 2013-2017 Martin Sandve Alnæs
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
from ffc.uflacs.integralgenerator import IntegralGenerator
from ffc.uflacs.language.format_lines import format_indented_lines


def generate_integral_code(ir, prefix, parameters):
    "Generate code for integral from intermediate representation."

    info("Generating code from ffc.uflacs representation")

    # FIXME: Is this the right precision value to use? Make it default to None or 0.
    precision = ir["integrals_metadata"]["precision"]

    # Create FFC C++ backend
    backend = FFCBackend(ir, parameters)

    # Configure kernel generator
    ig = IntegralGenerator(ir, backend, precision)

    # Generate code ast for the tabulate_tensor body
    parts = ig.generate()

    # Format code as string
    body = format_indented_lines(parts.cs_format(precision), 1)

    # Generate generic ffc code snippets and add uflacs specific parts
    code = initialize_integral_code(ir, prefix, parameters)
    code["tabulate_tensor"] = body
    code["additional_includes_set"] = set(ig.get_includes())
    code["additional_includes_set"].update(ir.get("additional_includes_set", ()))

    # TODO: Move to initialize_integral_code, this is not representation specific
    if ir.get("num_cells") is not None:
        ret = backend.language.Return(ir["num_cells"])
        code["num_cells"] = format_indented_lines(ret.cs_format(), 1)

    return code
