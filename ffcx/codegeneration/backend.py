# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Collection of FFCx specific pieces for the code generation phase."""

from __future__ import annotations

from ffcx.codegeneration.access import FFCXBackendAccess
from ffcx.codegeneration.definitions import FFCXBackendDefinitions
from ffcx.codegeneration.symbols import FFCXBackendSymbols
from ffcx.ir.representation import ExpressionIR, IntegralIR


class FFCXBackend:
    """Class collecting all aspects of the FFCx backend."""

    def __init__(self, ir: IntegralIR | ExpressionIR, options):
        """Initialise."""
        coefficient_numbering = ir.expression.coefficient_numbering
        coefficient_offsets = ir.expression.coefficient_offsets

        original_constant_offsets = ir.expression.original_constant_offsets

        self.symbols = FFCXBackendSymbols(
            coefficient_numbering, coefficient_offsets, original_constant_offsets
        )
        self.access = FFCXBackendAccess(
            ir.expression.entity_type, ir.expression.integral_type, self.symbols, options
        )
        self.definitions = FFCXBackendDefinitions(
            ir.expression.entity_type, ir.expression.integral_type, self.access, options
        )
