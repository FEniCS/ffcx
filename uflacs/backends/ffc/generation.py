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

from uflacs.generation.integralgenerator import IntegralGenerator

import uflacs.language.cnodes
from uflacs.language.format_lines import format_indented_lines
from uflacs.language.ufl_to_cnodes import UFL2CNodesTranslator
from uflacs.backends.ffc.access import FFCAccessBackend
from uflacs.backends.ffc.definitions import FFCDefinitionsBackend

class FFCBackend(object):
    "Class collecting all aspects of the FFC backend."
    def __init__(self, ir, parameters):
        self.language = uflacs.language.cnodes
        self.ufl_to_language = UFL2CNodesTranslator(self.language)
        self.definitions = FFCDefinitionsBackend(ir, self.language, parameters)
        self.access = FFCAccessBackend(ir, self.language, parameters)

def generate_tabulate_tensor_code(ir, prefix, parameters):

    # Create FFC C++ backend
    backend = FFCBackend(ir, parameters)

    # Create code generator for integral body
    ig = IntegralGenerator(ir, backend)

    # Generate code ast for the tabulate_tensor body
    parts = ig.generate()

    # Format code AST as one string
    body = format_indented_lines(parts.cs_format(), 1)
    #import IPython; IPython.embed()

    # Fetch includes
    includes = set()
    includes.update(ig.get_includes())
    includes.update(backend.definitions.get_includes())

    # Format uflacs specific code structures into a single
    # string and place in dict before returning to ffc
    code = {
        "tabulate_tensor": body,
        "additional_includes_set": includes,
    }
    return code
