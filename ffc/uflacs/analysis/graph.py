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

"""Linearized data structure for the computational graph."""

from ffc.uflacs.analysis.graph_vertices import build_graph_vertices
from ffc.uflacs.analysis.graph_symbols import build_graph_symbols


class Graph2(object):

    def __init__(self):

        self.nv = 0
        self.V = []
        self.e2i = {}

        self.expression_vertices = []

        self.V_shapes = []
        self.V_symbols = None  # Crs matrix
        self.total_unique_symbols = 0


def build_graph(expressions, DEBUG=False):

    # Make empty graph
    G = Graph2()

    # Populate with vertices
    G.e2i, G.V, G.expression_vertices = build_graph_vertices(expressions)
    G.nv = len(G.V)

    # Populate with symbols
    G.V_shapes, G.V_symbols, G.total_unique_symbols = \
        build_graph_symbols(G.V, G.e2i, DEBUG)

    if DEBUG:
        assert G.total_unique_symbols == len(set(G.V_symbols.data))

    return G
