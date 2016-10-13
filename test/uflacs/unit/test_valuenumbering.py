# -*- coding: utf-8 -*-

import pytest
from ufl import *
from ffc.uflacs.analysis.indexing import map_indexed_arg_components
from ffc.uflacs.analysis.indexing import map_component_tensor_arg_components

def test_map_index_arg_components():
    x = SpatialCoordinate(triangle)
    i, j, k = indices(3)

    # Tensors of rank 1, 2, 3
    A1 = x
    A2 = outer(x, x)
    A3 = outer(A2, x)

    indexed = A1[i]
    assert map_indexed_arg_components(indexed) == [0, 1]

    # Rank 2
    indexed = A2[i, j]
    assert map_indexed_arg_components(indexed) == [0, 1, 2, 3]
    indexed = A2[j, i]
    assert map_indexed_arg_components(indexed) == [0, 2, 1, 3]

    # Rank 3
    indexed = A3[i, j, k]
    assert map_indexed_arg_components(indexed) == [0, 1, 2, 3,
                                                   4, 5, 6, 7]
    indexed = A3[j, i, k]
    assert map_indexed_arg_components(indexed) == [0, 1, 4, 5,
                                                   2, 3, 6, 7]
    indexed = A3[i, k, j]
    assert map_indexed_arg_components(indexed) == [0, 2, 1, 3,
                                                   4, 6, 5, 7]
    indexed = A3[j, k, i]
    assert map_indexed_arg_components(indexed) == [0, 2, 4, 6,
                                                   1, 3, 5, 7]
    indexed = A3[k, i, j]
    assert map_indexed_arg_components(indexed) == [0, 4, 1, 5,
                                                   2, 6, 3, 7]
    indexed = A3[k, j, i]
    assert map_indexed_arg_components(indexed) == [0, 4, 2, 6,
                                                   1, 5, 3, 7]
    # Explanation:
    assert [ii*4+jj*2+kk  # strides are determined by relative position of i,j,k
            for kk in range(2)  # loop order matches indexing order
            for jj in range(2)  # loop range matches index dimensions
            for ii in range(2)] == [0, 4, 2, 6, 1, 5, 3, 7]


def test_map_component_tensor_arg_components():
    x = SpatialCoordinate(triangle)
    i, j, k = indices(3)

    # Tensors of rank 1, 2, 3
    A1 = x
    A2 = outer(x, x)
    A3 = outer(A2, x)
    Aij = A2[i,j]

    # Rank 1
    assert map_component_tensor_arg_components((x[i]*2.0)^(i,)) == [0,1]

    # Rank 2
    assert map_component_tensor_arg_components((A2[i,j]*2.0)^(i,j)) == [0,1,2,3]
    assert map_component_tensor_arg_components((A2[i,j]*2.0)^(j,i)) == [0,2,1,3]
    assert map_component_tensor_arg_components(A2[i,j]^(j,i)) == [0,2,1,3]

    # Rank 3
    assert map_component_tensor_arg_components((A3[i,j,k]*2.0)^(i,j,k)) == [0,1,2,3,4,5,6,7]
    assert map_component_tensor_arg_components((A3[i,j,k]*2.0)^(i,k,j)) == [0,2,1,3,4,6,5,7]
    assert map_component_tensor_arg_components((A3[i,j,k]*2.0)^(k,i,j)) == [0,2,4,6,1,3,5,7]

    # Explanation:
    assert [ii*4+jj*2+kk  # strides are determined by relative order of i,j,k in indexed expr
            for kk in range(2)  # loop order matches un-indexing order in component tensor
            for ii in range(2)  # loop range matches index dimensions
            for jj in range(2)] == [0,2,4,6,1,3,5,7]
