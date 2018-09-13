# -*- coding: utf-8 -*-
# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Assigning symbols to computational graph nodes."""

import numpy
from ufl import product

from ffc.uflacs.analysis.valuenumbering import ValueNumberer


def build_graph_symbols(V):
    """Tabulate scalar value numbering of all nodes in a a list based
    representation of an expression graph.

    Returns
    -------
    V_shapes - list of the total shapes of nodes in V.
    V_symbols - list of symbols (value numbers) of each component of each node in V.
    total_unique_symbols - The number of symbol values assigned to unique scalar
                           components of the nodes in V.

    """
    # Compute the total shape (value shape x index dimensions) for each node
    V_shapes = [(v.ufl_shape + v.ufl_index_dimensions) for v in V]

    # Compute the total value size for each node
    V_sizes = [product(sh) for sh in V_shapes]

    V_symbols = []
    value_numberer = ValueNumberer(V, V_sizes, V_symbols)

    for (i, v) in enumerate(V):
        V_symbols.append(value_numberer(v, i))

    return V_shapes, V_symbols, value_numberer.symbol_count
