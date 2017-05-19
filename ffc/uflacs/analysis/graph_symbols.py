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

"""Assigning symbols to computational graph nodes."""

import numpy
from ufl import product

from ffc.uflacs.analysis.crsarray import CRSArray
from ffc.uflacs.analysis.valuenumbering import ValueNumberer
from ffc.uflacs.analysis.expr_shapes import total_shape


def build_node_shapes(V):
    """Build total shapes for each node in list representation of expression graph.

    V is an array of ufl expressions, possibly nonscalar and with free indices.

    Returning a CRSArray where row i is the total shape of V[i].
    """
    # Dimensions of returned CRSArray
    nv = len(V)
    k = 0

    # Store shapes intermediately in an array of tuples
    V_shapes = numpy.empty(nv, dtype=object)
    for i, v in enumerate(V):
        # Compute total shape of V[i]
        tsh = total_shape(v)
        V_shapes[i] = tsh
        # Count number of elements for CRSArray representation
        k += len(tsh)

    # Return a more memory efficient CRSArray representation
    return CRSArray.from_rows(V_shapes, nv, k, int)


def build_node_sizes(V_shapes):
    "Compute all the products of a sequence of shapes."
    nv = len(V_shapes)
    V_sizes = numpy.zeros(nv, dtype=int)
    for i, sh in enumerate(V_shapes):
        V_sizes[i] = product(sh)
    return V_sizes


def build_node_symbols(V, e2i, V_shapes, V_sizes):
    """Tabulate scalar value numbering of all nodes in a a list based representation of an expression graph.

    Returns:
    V_symbols - CRSArray of symbols (value numbers) of each component of each node in V.
    total_unique_symbols - The number of symbol values assigned to unique scalar components of the nodes in V.
    """
    # "Sparse" int matrix for storing variable number of entries (symbols) per row (vertex),
    # with a capasity bounded by the number of scalar subexpressions including repetitions
    V_symbols = CRSArray(len(V), sum(V_sizes), int)

    # Visit each node with value numberer algorithm, storing the result for each as a row in the V_symbols CRSArray
    value_numberer = ValueNumberer(e2i, V_sizes, V_symbols)
    for i, v in enumerate(V):
        V_symbols.push_row(value_numberer(v, i))

    # Get the actual number of symbols created
    total_unique_symbols = value_numberer.symbol_count

    # assert all(x < total_unique_symbols for x in V_symbols.data)
    # assert (total_unique_symbols-1) in V_symbols.data

    return V_symbols, total_unique_symbols


def build_graph_symbols(V, e2i, DEBUG):
    """Tabulate scalar value numbering of all nodes in a a list based representation of an expression graph.

    Returns:
    V_shapes - CRSArray of the total shapes of nodes in V.
    V_symbols - CRSArray of symbols (value numbers) of each component of each node in V.
    total_unique_symbols - The number of symbol values assigned to unique scalar components of the nodes in V.
    """
    # Compute the total shape (value shape x index dimensions) for each node
    V_shapes = build_node_shapes(V)

    # Compute the total value size for each node
    V_sizes = build_node_sizes(V_shapes)

    # Mark values with symbols
    V_symbols, total_unique_symbols = build_node_symbols(V, e2i, V_shapes, V_sizes)

    return V_shapes, V_symbols, total_unique_symbols
