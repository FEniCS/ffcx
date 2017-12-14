# -*- coding: utf-8 -*-
"""
Tests of graph representation of expressions.
"""

from ufl import *
from ufl import product
from ufl.permutation import compute_indices

from ffc.uflacs.analysis.graph import build_graph
from ffc.uflacs.analysis.graph_rebuild import rebuild_expression_from_graph
# from ffc.uflacs.analysis.graph_rebuild import rebuild_scalar_e2i
# from ffc.uflacs.analysis.dependencies import (compute_dependencies,
#                                                mark_active,
#                                                mark_image)
# from ffc.uflacs.analysis.graph_ssa import (mark_partitions,
#                                       compute_dependency_count,
#                                       invert_dependencies,
#                                       default_cache_score_policy,
#                                       compute_cache_scores,
#                                       allocate_registers)

from operator import eq as equal


def test_graph_algorithm_allocates_correct_number_of_symbols():
    U = FiniteElement("CG", triangle, 1)
    V = VectorElement("CG", triangle, 1)
    W = TensorElement("CG", triangle, 1)
    u = Coefficient(U)
    v = Coefficient(V)
    w = Coefficient(W)

    # Testing some scalar expressions
    expr = u
    G = build_graph([expr], DEBUG=0)
    assert G.V_symbols.num_elements == 1
    assert G.total_unique_symbols == 1

    expr = u**2
    G = build_graph([expr], DEBUG=0)
    assert G.V_symbols.num_elements == 3
    assert G.total_unique_symbols == 3

    expr = u**2 / 2
    G = build_graph([expr], DEBUG=0)
    assert G.V_symbols.num_elements == 4
    assert G.total_unique_symbols == 4

    expr = dot(v, v) / 2
    G = build_graph([expr], DEBUG=0)
    assert G.V_symbols.num_elements == 5
    assert G.total_unique_symbols == 5

    # Testing Indexed
    expr = v[i] * v[i]
    G = build_graph([expr], DEBUG=0)
    assert G.V_symbols.num_elements == 2 + 2 + 2 + 1
    assert G.total_unique_symbols == 2 + 2 + 1

    # Reusing symbols for indexed with different ordering
    # Note that two index sums are created, giving 2+1 symbols
    expr = w[i, j] * w[j, i]
    G = build_graph([expr], DEBUG=0)
    assert G.V_symbols.num_elements == 4 + 4 + 4 + 4 + 2 + 1
    assert G.total_unique_symbols == 4 + 4 + 2 + 1

    # Testing ComponentTensor
    expr = dot(as_vector(2 * v[i], i), v)
    G = build_graph([expr], DEBUG=0)
    assert G.V_symbols.num_elements == 2 + 1 + 2 + 2 + 2 + 1
    assert G.total_unique_symbols == 2 + 1 + 2 + 1

    expr = dot(v + 2 * v, v)
    G = build_graph([expr], DEBUG=0)
    assert G.V_symbols.num_elements == 2 + 1 + 2 + 2 + 2 + 2 + 1
    assert G.total_unique_symbols == 2 + 1 + 2 + 2 + 1

    expr = outer(v, v)[i, j] * outer(v, v)[j, i]
    G = build_graph([expr], DEBUG=0)
    assert G.V_symbols.num_elements == 21  # 2+4+4+4 + 4+2+1
    assert G.total_unique_symbols == 13  # 2+4+4 + 2+1

    # Testing tensor/scalar
    expr = as_ufl(2)
    G = build_graph([expr], DEBUG=0)
    assert G.V_symbols.num_elements == 1
    assert G.total_unique_symbols == 1

    expr = v
    G = build_graph([expr], DEBUG=0)
    assert G.V_symbols.num_elements == 2
    assert G.total_unique_symbols == 2

    expr = outer(v, v)
    G = build_graph([expr], DEBUG=0)
    assert G.V_symbols.num_elements == 2 + 4
    assert G.total_unique_symbols == 2 + 4

    expr = as_tensor(v[i] * v[j], (i, j))
    G = build_graph([expr], DEBUG=0)
    assert G.V_symbols.num_elements == 2 + 2 + 2 + 4 + 4
    assert G.total_unique_symbols == 2 + 4

    expr = as_tensor(v[i] * v[j] / 2, (i, j))
    G = build_graph([expr], DEBUG=0)
    assert G.V_symbols.num_elements == 2 + 2 + 2 + 4 + 4 + 4 + 1
    assert G.total_unique_symbols == 2 + 1 + 4 + 4

    expr = outer(v, v) / 2  # converted to the tensor notation above
    G = build_graph([expr], DEBUG=0)
    assert G.V_symbols.num_elements == 2 + 2 + 2 + 4 + 4 + 4 + 1
    assert G.total_unique_symbols == 2 + 1 + 4 + 4


def test_rebuild_expression_from_graph_basic_scalar_expressions():
    U = FiniteElement("CG", triangle, 1)
    V = VectorElement("CG", triangle, 1)
    W = TensorElement("CG", triangle, 1)
    u = Coefficient(U)
    v = Coefficient(V)
    w = Coefficient(W)

    # ... Literals are reproduced
    literals = [as_ufl(0), as_ufl(1), as_ufl(3.14)]
    for v1 in literals:
        G = build_graph([v1])
        v2 = rebuild_expression_from_graph(G)
        assert v1 == v2

    v1 = u
    G = build_graph([v1])
    v2 = rebuild_expression_from_graph(G)
    assert v1 == v2

    # ... Simple operators are reproduced
    for v1 in [2 + u, u + u, u * u, u * 2, u**2, u**u, u / u, sin(u)]:
        G = build_graph([v1])
        v2 = rebuild_expression_from_graph(G)
        assert v1 == v2


def test_rebuild_expression_from_graph_on_products_with_indices():
    U = FiniteElement("CG", triangle, 1)
    V = VectorElement("CG", triangle, 1)
    W = TensorElement("CG", triangle, 1)
    u = Coefficient(U)
    v = Coefficient(V)
    w = Coefficient(W)
    i, j = indices(2)

    # Test fixed index
    fixed = [u * v[0], v[1] * v[0], w[0, 1] * w[0, 0]]
    for v1 in fixed:
        G = build_graph([v1])
        v2 = rebuild_expression_from_graph(G)
        assert v1 == v2

    # Test simple repeated index
    v1 = v[i] * v[i]
    G = build_graph([v1])
    v2 = rebuild_expression_from_graph(G)
    ve = v[0] * v[0] + v[1] * v[1]
    assert ve == v2

    # Test double repeated index
    v1 = w[i, j] * w[j, i]
    G = build_graph([v1])
    v2 = rebuild_expression_from_graph(G)
    ve = (w[1, 1] * w[1, 1] + w[1, 0] * w[0, 1]) + (w[0, 1] * w[1, 0] + w[0, 0] * w[0, 0])
    if 0:
        print()
        print(v1)
        print(ve)
        print(v2)
        print()
    assert ve == v2

    # Test mix of repeated and non-repeated index
    v1 = (w[i, j] * w[j, 0] + v[i]) * v[i]
    G = build_graph([v1])
    v2 = rebuild_expression_from_graph(G)
    ve = ((w[0, 0] * w[0, 0] + w[0, 1] * w[1, 0] + v[0]) * v[0]
          + (w[1, 0] * w[0, 0] + w[1, 1] * w[1, 0] + v[1]) * v[1])
    assert ve == v2


def test_rebuild_expression_from_graph_basic_tensor_expressions():
    U = FiniteElement("CG", triangle, 1)
    V = VectorElement("CG", triangle, 1)
    W = TensorElement("CG", triangle, 1)
    u = Coefficient(U)
    v = Coefficient(V)
    vb = Coefficient(V)
    w = Coefficient(W)
    wb = Coefficient(W)

    # Single vector
    v1 = v
    G = build_graph([v1])
    v2 = rebuild_expression_from_graph(G)
    assert as_vector((v1[0], v1[1])) == v2

    # Single tensor
    v1 = w
    G = build_graph([v1])
    v2 = rebuild_expression_from_graph(G)
    assert as_vector((v1[0, 0], v1[0, 1], v1[1, 0], v1[1, 1])) == v2

    # Vector sum
    v1 = v + v
    G = build_graph([v1])
    v2 = rebuild_expression_from_graph(G)
    ve = as_vector((v[0] + v[0], v[1] + v[1]))
    assert ve == v2

    v1 = v + vb
    G = build_graph([v1])
    v2 = rebuild_expression_from_graph(G)
    ve = as_vector((v[0] + vb[0], v[1] + vb[1]))
    assert ve == v2

    # Tensor sum
    v1 = w + w
    G = build_graph([v1])
    v2 = rebuild_expression_from_graph(G)
    ve = as_vector((w[0, 0] + w[0, 0], w[0, 1] + w[0, 1], w[1, 0] + w[1, 0], w[1, 1] + w[1, 1]))
    assert ve == v2

    v1 = w + wb
    G = build_graph([v1])
    v2 = rebuild_expression_from_graph(G)
    ve = as_vector((w[0, 0] + wb[0, 0], w[0, 1] + wb[0, 1], w[1, 0] + wb[1, 0], w[1, 1] + wb[1, 1]))
    assert ve == v2

    # Scalar-vector product
    v1 = u * v
    G = build_graph([v1])
    v2 = rebuild_expression_from_graph(G)
    ve = as_vector((u * v[0], u * v[1]))
    assert ve == v2

    # Scalar-tensor product
    v1 = u * w
    G = build_graph([v1])
    v2 = rebuild_expression_from_graph(G)
    ve = as_vector((u * w[0, 0], u * w[0, 1], u * w[1, 0], u * w[1, 1]))
    assert ve == v2

    # Vector-vector index based inner product
    v1 = v[i] * v[i]
    G = build_graph([v1])
    v2 = rebuild_expression_from_graph(G)
    ve = v[0] * v[0] + v[1] * v[1]
    assert ve == v2

    v1 = v[i] * vb[i]
    G = build_graph([v1])
    v2 = rebuild_expression_from_graph(G)
    ve = v[0] * vb[0] + v[1] * vb[1]
    assert ve == v2

    # Tensor-tensor index based transposed inner product
    v1 = w[i, j] * w[j, i]
    G = build_graph([v1])
    v2 = rebuild_expression_from_graph(G)
    ve = (w[0, 0] * w[0, 0] + w[0, 1] * w[1, 0]) \
        + (w[1, 0] * w[0, 1] + w[1, 1] * w[1, 1])
    assert ve == v2

    v1 = w[i, j] * wb[j, i]
    G = build_graph([v1])
    v2 = rebuild_expression_from_graph(G)
    ve = (w[0, 0] * wb[0, 0] + w[1, 0] * wb[0, 1]) \
        + (w[0, 1] * wb[1, 0] + w[1, 1] * wb[1, 1])
    assert ve == v2

    # Vector/scalar division
    v1 = v / u
    G = build_graph([v1])
    v2 = rebuild_expression_from_graph(G)
    ve = as_vector((v[0] / u, v[1] / u))
    assert ve == v2

    # Tensor/scalar division
    v1 = w / u
    G = build_graph([v1])
    v2 = rebuild_expression_from_graph(G)
    ve = as_vector((w[0, 0] / u, w[0, 1] / u, w[1, 0] / u, w[1, 1] / u))
    assert ve == v2

    # FIXME: Write more tests to discover bugs in ReconstructScalarSubexpressions.element_wise*

    # assert False

# Compounds not implemented, not expecting to do this anytime soon


def xtest_rebuild_expression_from_graph_on_compounds():
    U = FiniteElement("CG", triangle, 1)
    V = VectorElement("CG", triangle, 1)
    W = TensorElement("CG", triangle, 1)
    u = Coefficient(U)
    v = Coefficient(V)
    w = Coefficient(W)

    v1 = dot(v, v)
    G = build_graph([v1])
    v2 = rebuild_expression_from_graph(G)

    v1 = outer(v, v)
    G = build_graph([v1])
    v2 = rebuild_expression_from_graph(G)

    v1 = outer(v, v)[i, j] * outer(v, v)[j, i]
    G = build_graph([v1])
    v2 = rebuild_expression_from_graph(G)
    # print v1
    # print v2
    # FIXME: Assert something


def test_flattening_of_tensor_valued_expression_symbols():
    # from ffc.uflacs.analysis.graph import foo
    def flatten_expression_symbols(v, vsyms, opsyms):
        sh = v.ufl_shape
        if sh == ():
            assert len(vsyms) == 1
            if not opsyms:
                res = [(v, vsyms[0], ())]
            else:
                assert len(opsyms[0]) == 1
                res = [(v, vsyms[0], tuple(syms[0] for syms in opsyms))]
        else:
            res = []
            if isinstance(v, ufl.classes.Sum):
                for i in range(len(vsyms)):
                    u = None  # sum of component i for syms in opsyms
                    res += (u, vsyms[i], tuple(syms[i] for syms in opsyms))
        return res

    v = as_ufl(1)
    vsyms = (0,)
    opsyms = ()
    res = flatten_expression_symbols(v, vsyms, opsyms)
    assert res == [(v, 0, ())]
