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
    fi = v.free_indices()
    if fi:
        idims = v.index_dimensions()
        return tuple(idims[idx] for idx in sorted_by_count(fi))
    else:
        return ()


def compute_all_shapes(v):
    """Compute the tensor-, index-, and total shape of an expr.

    Returns (shape, size, index_shape, index_size, total_shape, total_size).
    """
    if isinstance(v, Condition):
        # TODO: Shape and index calls are invalid for conditions.
        #       Is this the best fix? Could also just return () from Condition.shape()?
        # Return scalar shape (conditions are scalar bool expressions so this works out fine)
        shape = ()
        index_shape = ()
        total_shape = ()
    else:
        shape = v.shape()
        index_shape = compute_index_shape(v)
        total_shape = shape + index_shape

    return (shape, index_shape, total_shape)


def total_shape(v):
    """Compute the total shape of an expr."""
    sh, ish, tsh = compute_all_shapes(v)
    return tsh
