# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Collection of FFCx specific pieces for the code generation phase."""

from ffcx.codegeneration.symbols import ReprManagerSymbols
import ffcx.codegeneration.lnodes as L
from ffcx.ir.representationutils import QuadratureRule
from ufl.core.expr import Expr


class ReprManager(object):
    """Class collecting all aspects of the FFCx backend."""

    def __init__(self, ir, options):

        coefficient_numbering = ir.coefficient_numbering
        coefficient_offsets = ir.coefficient_offsets

        original_constant_offsets = ir.original_constant_offsets

        self.ir = ir
        self.symbols = ReprManagerSymbols(coefficient_numbering,
                                          coefficient_offsets, original_constant_offsets)

        # TODO: Move this to symbols ?
        self.scopes = {}
        self.scopes = {quadrature_rule: {} for quadrature_rule in ir.integrand.keys()}
        self.scopes[None] = {}
        self.temp_symbols = {}

        # store quadrature rules indices
        self.quadrature_indices = {}
        self.quadrature_indices = {qr: self.create_quadrature_index(qr) for qr in ir.integrand.keys()}
        self.quadrature_indices[None] = self.create_quadrature_index(None)

    def get_var(self, quadrature_rule: QuadratureRule, v: Expr):
        """Lookup ufl expression v in variable scope dicts.

        Scope is determined by quadrature rule which identifies the
        quadrature loop scope or None if outside quadrature loops.

        If v is not found in quadrature loop scope, the piecewise
        scope (None) is checked.

        Returns the LNodes expression to access the value in the code.
        """
        assert isinstance(v, Expr)

        if v._ufl_is_literal_:
            return L.ufl_to_lnodes(v)

        # quadrature loop scope
        f = self.scopes[quadrature_rule].get(v)

        # piecewise scope
        if f is None:
            f = self.scopes[None].get(v)
        return f

    def set_var(self, quadrature_rule: QuadratureRule, v: Expr, vaccess: L.Symbol):
        """Set a new variable in variable scope dicts.

        Scope is determined by quadrature_rule which identifies the
        quadrature loop scope or None if outside quadrature loops.

        v is the ufl expression and vaccess is the LNodes
        expression to access the value in the code.

        """
        self.scopes[quadrature_rule][v] = vaccess

    def create_quadrature_index(self, quadrature_rule: QuadratureRule):
        """Return the quadrature index for the given quadrature rule."""

        ranges = [0]
        iq = self.symbols.quadrature_loop_index
        indices = [L.Symbol(iq.name, dtype=L.DataType.INT)]

        if quadrature_rule:
            ranges = [quadrature_rule.weights.size]
            if quadrature_rule.has_tensor_factors:
                dim = len(quadrature_rule.tensor_factors)
                ranges = [factor[1].size for factor in quadrature_rule.tensor_factors]
                indices = [L.Symbol(iq.name + f"{i}", dtype=L.DataType.INT) for i in range(dim)]

        return L.MultiIndex(indices, ranges)

    def create_dof_index(self, tabledata, index):
        """Create a multi index for the coefficient dofs."""
        name = index.name
        if tabledata.has_tensor_factorisation:
            dim = len(tabledata.tensor_factors)
            ranges = [factor.values.shape[-1] for factor in tabledata.tensor_factors]
            indices = [L.Symbol(f"{name}{i}", dtype=L.DataType.INT) for i in range(dim)]
        else:
            ranges = [tabledata.values.shape[-1]]
            indices = [L.Symbol(name, dtype=L.DataType.INT)]

        return L.MultiIndex(indices, ranges)
