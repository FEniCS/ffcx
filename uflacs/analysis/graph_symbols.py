
from ufl.common import sorted_by_count
from ufl.classes import Condition

from ffc.log import error

from uflacs.datastructures.arrays import int_array, object_array
from uflacs.datastructures.crs import CRS, rows_to_crs
from uflacs.analysis.valuenumbering import ValueNumberer


def total_shape(v):
    """Compute the total shape of an expr.

    The total shape is the regular shape tuple plus the index shape tuple.
    The index shape tuple is the tuple of index dimensions of the free indices
    of the expression, sorted by the count of the free indices.

    The total shape of a tensor valued expression A and A[*indices(A.rank())]
    is therefore the same.
    """
    if isinstance(v, Condition):
        # TODO: Shape and index calls are invalid for conditions.
        #       Is this the best fix? Could also just return () from Condition.shape()?
        tsh = ()

    else:
        # Regular shape
        sh = v.shape()

        # Index "shape"
        fi = v.free_indices()
        if fi:
            # Just an attempt at optimization, not running this code for expressions without free indices
            idims = v.index_dimensions()
            ish = tuple(idims[idx] for idx in sorted_by_count(fi))

            # Store "total" shape
            tsh = sh + ish

        else:
            tsh = sh

    return tsh


def build_node_shapes(V):
    """Build total shapes for each node in list representation of expression graph.

    V is an array of ufl expressions, possibly nonscalar and with free indices.

    Returning a CRS where row i is the total shape of V[i].
    """
    # Dimensions of returned CRS
    nv = len(V)
    k = 0

    # Store shapes intermediately in an array of tuples
    V_shapes = object_array(nv)
    for i, v in enumerate(V):
        # Compute total shape of V[i]
        V_shapes[i] = total_shape(v)

        # Count number of elements for CRS representation
        k += len(tsh)

    # Return a more memory efficient CRS representation
    return rows_to_crs(V_shapes, nv, k, int)


def build_node_sizes(V_shapes):
    "Compute all the products of a sequence of shapes."
    nv = len(V_shapes)
    V_sizes = int_array(nv)
    for i, sh in enumerate(V_shapes):
        V_sizes[i] = product(sh)
    return V_sizes


def build_node_symbols(V, e2i, V_shapes):
    """Tabulate scalar value numbering of all nodes in a a list based representation of an expression graph.

    Returns:
    V_symbols - CRS of symbols (value numbers) of each component of each node in V.
    total_unique_symbols - The number of symbol values assigned to unique scalar components of the nodes in V.
    """
    # Compute the total value size for each node, this gives an upper bound on the number of symbols we need
    V_sizes = build_node_sizes(V_shapes)
    max_symbols = sum(V_sizes)

    # "Sparse" int matrix for storing variable number of entries (symbols) per row (vertex).
    symbol_type = int
    V_symbols = CRS(len(V), max_symbols, symbol_type)

    # Visit each node with value numberer algorithm, storing the result for each as a row in the V_symbols CRS
    value_numberer = ValueNumberer(e2i, V_sizes, V_symbols)
    for i, v in enumerate(V):
        symbols = value_numberer.visit(v, i)
        V_symbols.push_row(symbols)

    total_unique_symbols = value_numberer.symbol_count

    assert all(x < total_unique_symbols for x in V_symbols.data)
    assert (total_unique_symbols-1) in V_symbols.data

    return V_symbols, total_unique_symbols


def build_graph_symbols(V, e2i, DEBUG):
    """Tabulate scalar value numbering of all nodes in a a list based representation of an expression graph.

    Returns:
    V_shapes - CRS of the total shapes of nodes in V.
    V_symbols - CRS of symbols (value numbers) of each component of each node in V.
    total_unique_symbols - The number of symbol values assigned to unique scalar components of the nodes in V.
    """
    # Compute the total shape (value shape x index dimensions) for each node
    V_shapes = build_node_shapes(V)

    # Mark values with symbols
    V_symbols, total_unique_symbols = build_node_symbols(V, e2i, V_shapes)

    return V_shapes, V_symbols, total_unique_symbols
