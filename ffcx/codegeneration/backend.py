# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Collection of FFCx specific pieces for the code generation phase."""

from ffcx.codegeneration.access import FFCXBackendAccess
from ffcx.codegeneration.definitions import FFCXBackendDefinitions
from ffcx.codegeneration.symbols import FFCXBackendSymbols
import ffcx.codegeneration.lnodes as L


class FFCXBackend(object):
    """Class collecting all aspects of the FFCx backend."""

    def __init__(self, ir, options):

        coefficient_numbering = ir.coefficient_numbering
        coefficient_offsets = ir.coefficient_offsets

        original_constant_offsets = ir.original_constant_offsets

        self.symbols = FFCXBackendSymbols(coefficient_numbering,
                                          coefficient_offsets, original_constant_offsets)
        self.access = FFCXBackendAccess(ir, self.symbols, options)
        self.definitions = FFCXBackendDefinitions(ir, self.access, options)

        self.scopes = {}
        self.scopes = {quadrature_rule: {} for quadrature_rule in ir.integrand.keys()}
        self.scopes[None] = {}

    def get_var(self, quadrature_rule, v):
        """Lookup ufl expression v in variable scope dicts.

        Scope is determined by quadrature rule which identifies the
        quadrature loop scope or None if outside quadrature loops.

        If v is not found in quadrature loop scope, the piecewise
        scope (None) is checked.

        Returns the LNodes expression to access the value in the code.
        """
        if v._ufl_is_literal_:
            return L.ufl_to_lnodes(v)

        # quadrature loop scope
        f = self.scopes[quadrature_rule].get(v)

        # piecewise scope
        if f is None:
            f = self.scopes[None].get(v)
        return f

    def set_var(self, quadrature_rule, v, vaccess):
        """Set a new variable in variable scope dicts.

        Scope is determined by quadrature_rule which identifies the
        quadrature loop scope or None if outside quadrature loops.

        v is the ufl expression and vaccess is the LNodes
        expression to access the value in the code.

        """
        self.scopes[quadrature_rule][v] = vaccess
