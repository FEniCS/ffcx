# -*- coding: utf-8 -*-
# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Linearized data structure for the computational graph."""

from ffc.uflacs.analysis.graph_vertices import build_graph_vertices
from ffc.uflacs.analysis.graph_symbols import build_graph_symbols


class Graph2(object):
    def __init__(self):

        self.V = []
        self.e2i = {}

        self.expression_vertices = []

        self.V_shapes = []
        self.V_symbols = None
        self.total_unique_symbols = 0


def build_graph(expressions):

    G = Graph2()

    # Populate with vertices
    G.e2i, G.V, G.expression_vertices = build_graph_vertices(expressions, scalar=False)

    # Populate with symbols
    G.V_shapes, G.V_symbols, G.total_unique_symbols = build_graph_symbols(G.V)

    return G
