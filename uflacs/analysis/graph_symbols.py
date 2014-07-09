
from ufl.classes import (Terminal, FormArgument, Grad, Restricted,
                         Indexed, ComponentTensor, ListTensor, Transposed, Variable,
                         IndexSum, MultiIndex, Condition,
                         UtilityType, Label, ExprList, ExprMapping)

from uflacs.utils.log import error

from uflacs.analysis.datastructures import (int_array, object_array,
                                              CRS, rows_to_crs, rows_dict_to_crs)

from uflacs.analysis.indexing import *


def build_node_shapes(V):
    nv = len(V)
    k = 0
    V_shapes = object_array(nv)
    for i,v in enumerate(V):

        if isinstance(v, Condition):
            # FIXME: Shape and index calls are invalid for conditions. Is this the best fix?
            tsh = ()
        else:
            # Regular shape
            sh = v.shape()
            # Index "shape"
            idims = v.index_dimensions()
            ish = tuple(idims[idx] for idx in sorted_indices(v.free_indices()))
            # Store "total" shape and size
            tsh = sh + ish

        V_shapes[i] = tsh
        # Count number of elements for CRS representation
        k += len(tsh)

    # Return a more memory efficient CRS representation
    return rows_to_crs(V_shapes, nv, k, int)

def build_node_sizes(V_shapes):
    nv = len(V_shapes)
    V_sizes = int_array(nv)
    for i,sh in enumerate(V_shapes):
        V_sizes[i] = product(sh)
    return V_sizes

def get_node_symbols(expr, e2i, V_symbols):
    return V_symbols[e2i[expr]]

def map_indexed_symbols(v, e2i, V_symbols):
    # Reuse symbols of arg A for Aii
    Aii = v
    A = Aii.operands()[0]

    # Get symbols of argument A
    A_symbols = get_node_symbols(A, e2i, V_symbols)

    # Map A_symbols to Aii_symbols
    d = map_indexed_arg_components(Aii)
    symbols = [A_symbols[k] for k in d]
    return symbols

def map_component_tensor_symbols(v, e2i, V_symbols):
    # Reuse symbols of arg Aii for A
    A = v
    Aii = A.operands()[0]

    # Get symbols of argument Aii
    Aii_symbols = get_node_symbols(Aii, e2i, V_symbols)

    # Map A_symbols to Aii_symbols
    d = map_component_tensor_arg_components(A)
    symbols = [Aii_symbols[k] for k in d]
    return symbols

def map_list_tensor_symbols(v, e2i, V_symbols):
    A = v
    rows = A.operands()

    row_symbols = [get_node_symbols(row, e2i, V_symbols) for row in rows]

    symbols = []
    for rowsymb in row_symbols:
        symbols.extend(rowsymb) # FIXME: Test that this produces the right transposition
    return symbols

def map_transposed_symbols(v, e2i, V_symbols):
    AT = v
    A, = AT.operands()

    assert not A.free_indices(), "Assuming no free indices in transposed (for now), report as bug if needed." # FIXME
    r, c = A.shape()

    A_symbols = get_node_symbols(A, e2i, V_symbols)
    assert len(A_symbols) == r*c

    # AT[j,i] = A[i,j]
    # sh(A) = (r,c)
    # sh(AT) = (c,r)
    # AT[j*r+i] = A[i*c+j]
    symbols = [None]*(r*c)
    for j in xrange(c):
        for i in xrange(r):
            symbols[j*r+i] = A_symbols[i*c+j]

    return symbols

def map_variable_symbols(v, e2i, V_symbols):
    # Direct reuse of all symbols
    return get_node_symbols(v.operands()[0], e2i, V_symbols)

mappable_type = (Indexed, ComponentTensor, ListTensor, Transposed, Variable)
def map_symbols(v, e2i, V_symbols):
    if isinstance(v, Indexed):
        symbols = map_indexed_symbols(v, e2i, V_symbols)
    elif isinstance(v, ComponentTensor):
        symbols = map_component_tensor_symbols(v, e2i, V_symbols)
    elif isinstance(v, ListTensor):
        symbols = map_list_tensor_symbols(v, e2i, V_symbols)
    elif isinstance(v, Transposed):
        symbols = map_transposed_symbols(v, e2i, V_symbols)
    elif isinstance(v, Variable):
        symbols = map_variable_symbols(v, e2i, V_symbols)
    else:
        error("Not a mappable type!")
    return symbols

def build_node_symbols(V, e2i, V_shapes):

    # Compute the total value size for each node, this gives the max number of symbols we need
    V_sizes = build_node_sizes(V_shapes)
    max_symbols = sum(V_sizes)

    # Sparse int matrix for storing variable number of entries (symbols) per row (vertex).
    symbol_type = int
    V_symbols = CRS(len(V), max_symbols, symbol_type)

    # Generator for new symbols with a running counter
    def new_symbols(n):
        a = new_symbols.symbol_count
        b = a + n
        new_symbols.symbol_count = b
        return xrange(a, b)
    new_symbols.symbol_count = 0

    # For all vertices
    for i,v in enumerate(V):
        n = V_sizes[i]

        if isinstance(v, mappable_type):
            # Map symbols for expressions that only represent a different
            # view of other expressions through shape and indexing mappings.
            symbols = map_symbols(v, e2i, V_symbols)
        elif isinstance(v, FormArgument):
            # Create new symbols for expressions that represent new values
            # TODO: Ignoring symmetries for now, handle by creating only
            # some new symbols and mapping the rest using the symmetry map.
            symbols = new_symbols(n)
        else:
            # Create new symbols for expressions that represent new values
            symbols = new_symbols(n)

        assert len(symbols) == n
        V_symbols.push_row(symbols)

    total_unique_symbols = new_symbols.symbol_count

    assert all(x < total_unique_symbols for x in V_symbols.data)
    assert (total_unique_symbols-1) in V_symbols.data

    return V_symbols, total_unique_symbols

def build_graph_symbols(V, e2i, DEBUG):
    # Compute the total shape (value shape x index dimensions) for each node
    V_shapes = build_node_shapes(V)

    # Mark values with symbols
    V_symbols, total_unique_symbols = build_node_symbols(V, e2i, V_shapes)

    return V_shapes, V_symbols, total_unique_symbols
