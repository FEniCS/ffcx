__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-11-03"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# Python modules
import Numeric

def build_indices(dims):
    """Create a list of all index combinations matching the given list
    of index dimensions. Someone please tell me if there is a better
    way to iterate over a multi-dimensional array. The list of indices
    is constucted by iterating over all integer numbers that can be
    represented with the base of each position determined by the given
    dimension for that index."""
    if not dims:
        return []
    current = Numeric.zeros(len(dims))
    indices = []
    posvalue = [1] + list(Numeric.cumproduct(dims)[:-1])
    for i in range(Numeric.product(dims)):
        for pos in range(len(dims)):
            current[len(dims) - 1 - pos] = (i / posvalue[pos]) % dims[pos]
        indices += [list(current)]
    return indices

class MultiIndex:

    """A MultiIndex represents a list of indices and holds the
    following data:

        rank    - rank of multiindex
        dims    - a list of dimensions
        indices - list of all possible multiindex values"""

    def __init__(self, dims):
        "Create MultiIndex from given list of dimensions."
        self.rank = len(dims)
        self.dims = [] + dims
        self.indices = build_indices(dims)
        return

    def __repr__(self):
        "Print nicely formatted representation of MultiIndex."
        return "rank = %d dims = %s" % (self.rank, str(self.dims))
