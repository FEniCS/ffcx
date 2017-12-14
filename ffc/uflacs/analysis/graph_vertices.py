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
# along with UFLACS. If not, see <http://www.gnu.org/licenses/>.

"""Algorithms for working with graphs."""

import numpy

from ufl.classes import MultiIndex, Label

from ffc.uflacs.analysis.modified_terminals import is_modified_terminal


def count_nodes_with_unique_post_traversal(expr, e2i=None, skip_terminal_modifiers=False):
    """Yields o for each node o in expr, child before parent.
    Never visits a node twice."""
    if e2i is None:
        e2i = {}

    def getops(e):
        "Get a modifyable list of operands of e, optionally treating modified terminals as a unit."
        # TODO: Maybe use e._ufl_is_terminal_modifier_
        if e._ufl_is_terminal_ or (skip_terminal_modifiers and is_modified_terminal(e)):
            return []
        else:
            return list(e.ufl_operands)

    stack = []
    stack.append((expr, getops(expr)))
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


def build_array_from_counts(e2i):
    nv = len(e2i)
    V = numpy.empty(nv, dtype=object)
    for e, i in e2i.items():
        V[i] = e
    return V


def build_node_counts(expressions):
    e2i = {}
    for expr in expressions:
        count_nodes_with_unique_post_traversal(expr, e2i, False)
    return e2i


def build_scalar_node_counts(expressions):
    # Count unique expression nodes across multiple expressions
    e2i = {}
    for expr in expressions:
        count_nodes_with_unique_post_traversal(expr, e2i, True)
    return e2i


def build_graph_vertices(expressions):
    # Count unique expression nodes
    e2i = build_node_counts(expressions)

    # Make a list of the nodes by their ordering
    V = build_array_from_counts(e2i)

    # Get vertex indices representing input expression roots
    ri = [e2i[expr] for expr in expressions]

    return e2i, V, ri


def build_scalar_graph_vertices(expressions):
    # Count unique expression nodes across multiple expressions, treating modified terminals as a unit
    e2i = build_scalar_node_counts(expressions)

    # Make a list of the nodes by their ordering
    V = build_array_from_counts(e2i)

    # Get vertex indices representing input expression roots
    ri = [e2i[expr] for expr in expressions]

    return e2i, V, ri
