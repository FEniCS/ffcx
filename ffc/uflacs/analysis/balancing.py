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

"""Algorithms for the representation phase of the form compilation."""

from ufl.classes import (ReferenceValue, ReferenceGrad, Grad,
                         CellAvg, FacetAvg,
                         PositiveRestricted, NegativeRestricted,
                         Indexed)
from ufl.corealg.multifunction import MultiFunction
from ufl.corealg.map_dag import map_expr_dag


modifier_precedence = [ReferenceValue, ReferenceGrad, Grad,
                       CellAvg, FacetAvg,
                       PositiveRestricted, NegativeRestricted,
                       Indexed]

modifier_precedence = { m._ufl_handler_name_: i for i, m in enumerate(modifier_precedence) }

# TODO: Move this to ufl?
# TODO: Add expr._ufl_modifier_precedence_ ? Add Terminal first and Operator last in the above list.


def balance_modified_terminal(expr):
    # NB! Assuminge e.g. grad(cell_avg(expr)) does not occur,
    # i.e. it is simplified to 0 immediately.

    if expr._ufl_is_terminal_:
        return expr

    assert expr._ufl_is_terminal_modifier_

    orig = expr

    # Build list of modifier layers
    layers = [expr]
    while not expr._ufl_is_terminal_:
        assert expr._ufl_is_terminal_modifier_
        expr = expr.ufl_operands[0]
        layers.append(expr)
    assert layers[-1] is expr
    assert expr._ufl_is_terminal_

    # Apply modifiers in order
    layers = sorted(layers[:-1], key=lambda e: modifier_precedence[e._ufl_handler_name_])
    for op in layers:
        ops = (expr,) + op.ufl_operands[1:]
        expr = op._ufl_expr_reconstruct_(*ops)

    # Preserve id if nothing has changed
    return orig if expr == orig else expr


class BalanceModifiers(MultiFunction):

    def expr(self, expr, *ops):
        return expr._ufl_expr_reconstruct_(*ops)

    def terminal(self, expr):
        return expr

    def _modifier(self, expr, *ops):
        return balance_modified_terminal(expr)

    reference_value = _modifier
    reference_grad = _modifier
    grad = _modifier
    cell_avg = _modifier
    facet_avg = _modifier
    positive_restricted = _modifier
    negative_restricted = _modifier


def balance_modifiers(expr):
    mf = BalanceModifiers()
    return map_expr_dag(mf, expr)
