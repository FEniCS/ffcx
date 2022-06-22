# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Collection of FFCx specific pieces for the code generation phase."""

import ffcx.codegeneration.C.cnodes
from ffcx.codegeneration.access import FFCXBackendAccess
from ffcx.codegeneration.C.ufl_to_cnodes import UFL2CNodesTranslatorCpp
from ffcx.codegeneration.definitions import FFCXBackendDefinitions
from ffcx.codegeneration.symbols import FFCXBackendSymbols


class FFCXBackend(object):
    """Class collecting all aspects of the FFCx backend."""

    def __init__(self, ir, parameters):

        # This is the seam where cnodes/C is chosen for the FFCx backend
        self.language = ffcx.codegeneration.C.cnodes
        scalar_type = parameters["scalar_type"]
        self.ufl_to_language = UFL2CNodesTranslatorCpp(self.language, scalar_type)

        coefficient_numbering = ir.coefficient_numbering
        coefficient_offsets = ir.coefficient_offsets

        original_constant_offsets = ir.original_constant_offsets

        self.symbols = FFCXBackendSymbols(self.language, coefficient_numbering,
                                          coefficient_offsets, original_constant_offsets)
        self.definitions = FFCXBackendDefinitions(ir, self.language,
                                                  self.symbols, parameters)
        self.access = FFCXBackendAccess(ir, self.language, self.symbols,
                                        parameters)
