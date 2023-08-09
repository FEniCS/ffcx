# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Collection of FFCx specific pieces for the code generation phase."""

from ffcx.codegeneration.access import FFCXAccess
from ffcx.codegeneration.definitions import FFCXDefinitions
from ffcx.codegeneration.symbols import FFCXSymbols


class FFCXBackend(object):
    """Class collecting all aspects of the FFCx backend."""

    def __init__(self, ir, options):
        self.symbols = FFCXSymbols(ir)
        self.definitions = FFCXDefinitions(ir, self.symbols, options)
        self.access = FFCXAccess(ir, self.symbols, options)
