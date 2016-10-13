# -*- coding: utf-8 -*-
"""Tests of utilities for dealing with ufl indexing and components vs
flattened index spaces.

"""

from ufl import *
from ufl import product
from ufl.permutation import compute_indices

from ffc.uflacs.analysis.indexing import (map_indexed_arg_components,
                                      map_component_tensor_arg_components)
from ffc.uflacs.analysis.graph_symbols import (map_list_tensor_symbols,
                                           map_transposed_symbols,
                                           get_node_symbols)
from ffc.uflacs.analysis.graph import build_graph

from operator import eq as equal


def test_map_indexed_arg_components():
    W = TensorElement("CG", triangle, 1)
    A = Coefficient(W)
    i, j = indices(2)

    # Ordered indices:
    d = map_indexed_arg_components(A[i, j])
    assert equal(d, [0, 1, 2, 3])

    # Swapped ordering of indices:
    d = map_indexed_arg_components(A[j, i])
    assert equal(d, [0, 2, 1, 3])


def test_map_indexed_arg_components2():

    # This was the previous return type, copied here to preserve the
    # test without having to rewrite
    def map_indexed_arg_components2(Aii):
        c1, c2 = map_indexed_to_arg_components(Aii)
        d = [None] * len(c1)
        for k in range(len(c1)):
            d[c1[k]] = k
        return d

    W = TensorElement("CG", triangle, 1)
    A = Coefficient(W)
    i, j = indices(2)

    # Ordered indices:
    d = map_indexed_arg_components2(A[i, j])
    assert equal(d, [0, 1, 2, 3])

    # Swapped ordering of indices:
    d = map_indexed_arg_components2(A[j, i])
    assert equal(d, [0, 2, 1, 3])


def test_map_componenttensor_arg_components():
    W = TensorElement("CG", triangle, 1)
    A = Coefficient(W)
    i, j = indices(2)

    # Ordered indices:
    d = map_component_tensor_arg_components(as_tensor(2 * A[i, j], (i, j)))
    assert equal(d, [0, 1, 2, 3])

    # Swapped ordering of indices:
    d = map_component_tensor_arg_components(as_tensor(2 * A[i, j], (j, i)))
    assert equal(d, [0, 2, 1, 3])


def test_map_list_tensor_symbols():
    U = FiniteElement("CG", triangle, 1)
    u = Coefficient(U)
    A = as_tensor(((u + 1, u + 2, u + 3), (u**2 + 1, u**2 + 2, u**2 + 3)))
    # Would be nicer to refactor build_graph a bit so we could call
    # map_list_tensor_symbols directly...
    G = build_graph([A], DEBUG=False)
    s1 = list(get_node_symbols(A, G.e2i, G.V_symbols))
    s2 = [get_node_symbols(e, G.e2i, G.V_symbols)[0] for e in (u + 1, u + 2, u + 3, u**2 + 1, u**2 + 2, u**2 + 3)]
    assert s1 == s2


def test_map_transposed_symbols():
    W = TensorElement("CG", triangle, 1)
    w = Coefficient(W)
    A = w.T
    # Would be nicer to refactor build_graph a bit so we could call
    # map_transposed_symbols directly...
    G = build_graph([A], DEBUG=False)
    s1 = list(get_node_symbols(A, G.e2i, G.V_symbols))
    s2 = list(get_node_symbols(w, G.e2i, G.V_symbols))
    s2[1], s2[2] = s2[2], s2[1]
    assert s1 == s2

    W = TensorElement("CG", tetrahedron, 1)
    w = Coefficient(W)
    A = w.T
    # Would be nicer to refactor build_graph a bit so we could call
    # map_transposed_symbols directly...
    G = build_graph([A], DEBUG=False)
    s1 = list(get_node_symbols(A, G.e2i, G.V_symbols))
    s2 = list(get_node_symbols(w, G.e2i, G.V_symbols))
    s2[1], s2[2], s2[5], s2[3], s2[6], s2[7] = s2[3], s2[6], s2[7], s2[1], s2[2], s2[5]
    assert s1 == s2
