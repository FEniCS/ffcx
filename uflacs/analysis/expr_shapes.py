"""Tools for computing various shapes of ufl expressions.

The total shape is the regular shape tuple plus the index shape tuple.
The index shape tuple is the tuple of index dimensions of the free indices
of the expression, sorted by the count of the free indices.

The total shape of a tensor valued expression A and A[*indices(A.rank())]
is therefore the same.
"""

from ufl.common import sorted_by_count
from ufl.classes import Condition


def compute_index_shape(v):
    "Compute the 'index shape' of v."
    return v.ufl_index_dimensions


def compute_all_shapes(v):
    """Compute the tensor-, index-, and total shape of an expr.

    Returns (shape, size, index_shape, index_size, total_shape, total_size).
    """
    shape = v.ufl_shape
    index_shape = v.ufl_index_dimensions
    total_shape = shape + index_shape
    return (shape, index_shape, total_shape)


def total_shape(v):
    """Compute the total shape of an expr."""
    sh, ish, tsh = compute_all_shapes(v)
    return tsh
