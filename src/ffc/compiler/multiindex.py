__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-03 -- 2005-09-13"
__copyright__ = "Copyright (C) 2004, 2005 Anders Logg"
__license__  = "GNU GPL Version 2"

# Python modules
import Numeric

# FFC common modules
from ffc.common.util import *

def build_indices(dims):
    """Create a list of all index combinations matching the given list
    of index dimensions. Someone please tell me if there is a better
    way to iterate over a multi-dimensional array. The list of indices
    is constucted by iterating over all integer numbers that can be
    represented with the base of each position determined by the given
    dimension for that index."""
    if not dims:
        return []
    rdims = [] + dims;
    rdims.reverse()
    current = Numeric.zeros(len(rdims))
    indices = []
    posvalue = [1] + list(Numeric.cumproduct(rdims)[:-1])
    for i in range(Numeric.product(rdims)):
        for pos in range(len(rdims)):
            j = len(rdims) - 1
            current[len(rdims) - 1 - pos] = (i / posvalue[pos]) % rdims[pos]
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
        self.dims = listcopy(dims)
        self.indices = build_indices(dims)
        return

    def __repr__(self):
        "Print nicely formatted representation of MultiIndex."
        return "rank = %d dims = %s" % (self.rank, str(self.dims))
