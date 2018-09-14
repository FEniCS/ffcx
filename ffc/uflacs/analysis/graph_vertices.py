# -*- coding: utf-8 -*-
# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Algorithms for working with graphs."""

from ffc.uflacs.analysis.modified_terminals import is_modified_terminal
from ufl.classes import Label, MultiIndex


def count_nodes_with_unique_post_traversal(expr, e2i=None,
                                           skip_terminal_modifiers=False):
    """Yields o for each node o in expr, child before parent.
    Never visits a node twice."""
    if e2i is None:
        e2i = {}

    def getops(e):
        """Get a modifiable list of operands of e, optionally treating modified terminals as a unit."""
        # TODO: Maybe use e._ufl_is_terminal_modifier_
        if e._ufl_is_terminal_ or (skip_terminal_modifiers
                                   and is_modified_terminal(e)):
            return []
        else:
            return list(e.ufl_operands)

    stack = [(expr, getops(expr))]
    while stack:
        expr, ops = stack[-1]
        for i, o in enumerate(ops):
            if o is not None and o not in e2i:
                stack.append((o, getops(o)))
                ops[i] = None
                break
        else:
            if not isinstance(expr, (MultiIndex, Label)):
                count = len(e2i)
                e2i[expr] = count
            stack.pop()
    return e2i


def build_graph_vertices(expressions, scalar=False):
    # Count unique expression nodes

    e2i = {}
    for expr in expressions:
        count_nodes_with_unique_post_traversal(expr, e2i, scalar)

    # Invert the map to get index->expression
    V = sorted(e2i, key=e2i.get)

    # Get vertex indices representing input expression roots
    expression_vertices = [e2i[expr] for expr in expressions]

    return e2i, V, expression_vertices
