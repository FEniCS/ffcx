# -*- coding: utf-8 -*-
# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Collection of FFC specific pieces for the code generation phase."""

import ffc.codegeneration.C.cnodes
from ffc.codegeneration.access import FFCBackendAccess
from ffc.codegeneration.C.ufl_to_cnodes import UFL2CNodesTranslatorCpp
from ffc.codegeneration.definitions import FFCBackendDefinitions
from ffc.codegeneration.symbols import FFCBackendSymbols


class FFCBackend(object):
    """Class collecting all aspects of the FFC backend."""

    def __init__(self, ir, parameters):

        # This is the seam where cnodes/C is chosen for the ffc backend
        self.language = ffc.codegeneration.C.cnodes
        scalar_type = parameters.get("scalar_type", "double")
        self.ufl_to_language = UFL2CNodesTranslatorCpp(self.language, scalar_type)

        coefficient_numbering = ir.coefficient_numbering
        coefficient_offsets = ir.coefficient_offsets
        self.symbols = FFCBackendSymbols(self.language, coefficient_numbering,
                                         coefficient_offsets)
        self.definitions = FFCBackendDefinitions(ir, self.language,
                                                 self.symbols, parameters)
        self.access = FFCBackendAccess(ir, self.language, self.symbols,
                                       parameters)
