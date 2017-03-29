# -*- coding: utf-8 -*-
# Copyright (C) 2011-2017 Martin Sandve Alnæs
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

"""Collection of FFC specific pieces for the code generation phase."""

import ffc.uflacs.language.cnodes
from ffc.uflacs.language.ufl_to_cnodes import UFL2CNodesTranslatorCpp

from ffc.uflacs.backends.ffc.symbols import FFCBackendSymbols
from ffc.uflacs.backends.ffc.access import FFCBackendAccess
from ffc.uflacs.backends.ffc.definitions import FFCBackendDefinitions


class FFCBackend(object):
    "Class collecting all aspects of the FFC backend."
    def __init__(self, ir, parameters):

        # This is the seam where cnodes/C++ is chosen for the ffc backend
        self.language = ffc.uflacs.language.cnodes
        self.ufl_to_language = UFL2CNodesTranslatorCpp(self.language)

        coefficient_numbering = ir["coefficient_numbering"]
        self.symbols = FFCBackendSymbols(self.language, coefficient_numbering)
        self.definitions = FFCBackendDefinitions(ir, self.language, self.symbols, parameters)
        self.access = FFCBackendAccess(ir, self.language, self.symbols, parameters)
