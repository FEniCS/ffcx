__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-03"
__copyright__ = "Copyright (C) 2004-2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Garth N. Wells 2006
# Modified by Marie E. Rognes (meg@math.uio.no) 2007
# Modified by Kristian B. Oelgaard, 2009
# Last changed: 2009-12-09

# Python modules.
import numpy

# FFC modules.
from ffc.utils import listcopy
from ffc.log import error

def build_indices(dims):
    "Create a list of all index combinations"
    if not dims: return [[]]
    ranges = listcopy(dims)
    return reduce(outer_join, ranges, [[]])

def outer_join(a, b):
    """Let a be a list of lists and b a list. We append each element
    of b to each list in a and return the resulting list of lists."""
    outer = [] 
    for i in range(len(a)):
        for j in range(len(b)):
            outer += [a[i] + [b[j]]]
    return outer

def create_multi_index(monomial, index_type):
    "Find dimensions and create multi index for monomial of given index type."

    # Get sorted unique monomial indices
    indices = monomial.extract_unique_indices(index_type)
    indices.sort()

    # Check that we got all indices correctly
    for (i, index) in enumerate(indices):
        if not i == index.index_id:
            error("Unable to extract all indices.")

    # Get dimensions
    dims = [range(len(index.index_range)) for index in indices]

    return MultiIndex(dims)

class MultiIndex:

    """A MultiIndex represents a list of indices and holds the
    following data:

        rank    - rank of multiindex
        dims    - a list of dimensions
        indices - list of all possible multiindex values"""

    def __init__(self, dims):
        "Create multi index from given list of ranges"
        self.rank = len(dims)
        self.dims = [len(dim) for dim in dims]
        self.indices = build_indices(dims)
        return

    def __repr__(self):
        "Pretty print"
        return "rank = %d dims = %s indices = %s" % (self.rank, str(self.dims), str(self.indices))
