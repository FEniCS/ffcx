# -*- coding: utf-8 -*-
# Copyright (C) 2011-2015 Martin Sandve Alnæs
#
# This file is part of UFLACS.
#
# UFLACS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# UFLACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with UFLACS. If not, see <http://www.gnu.org/licenses/>

"""FFC specific algorithms for the generation phase."""

from uflacs.codeutils.cpp_expr_formatting_rules import CppExprFormatter
from uflacs.generation.integralgenerator import IntegralGenerator
from uflacs.backends.ffc.access import FFCAccessBackend
from uflacs.backends.ffc.definitions import FFCDefinitionsBackend


def generate_tabulate_tensor_code(ir, parameters):

    # Create C++ backend
    language_formatter = CppExprFormatter()

    # Create FFC backend
    backend_access = FFCAccessBackend(ir, parameters)
    backend_definitions = FFCDefinitionsBackend(ir, parameters)

    # Create code generator for integral body
    ig = IntegralGenerator(ir, language_formatter, backend_access, backend_definitions)

    # Generate code for the tabulate_tensor body
    body = ig.generate()

    # Fetch includes
    includes = set()
    includes.update(ig.get_includes())
    includes.update(backend_definitions.get_includes())

    # Format uflacs specific code structures into a single
    # string and place in dict before returning to ffc
    code = {
        "tabulate_tensor": body,
        "additional_includes_set": includes,
    }
    return code
