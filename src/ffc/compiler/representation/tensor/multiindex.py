__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-03 -- 2007-03-08"
__copyright__ = "Copyright (C) 2004-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# Modified by Garth N. Wells 2006
# Modified by Marie E. Rognes (meg@math.uio.no) 2007

# Python modules
import numpy

# FFC common modules
from ffc.common.utils import *

def build_indices(dims):
    "Create a list of all index combinations"
    if not dims:
        return [[]]
    dims.insert(0, [[]]) # Special first case.
    return reduce(outer_join, dims)

def outer_join(a, b):
    """Let a be a list of lists and b a list. We append each element
    of b to each list in a and return the resulting list of lists."""
    outer = [] 
    for i in range(len(a)):
        for j in range(len(b)):
            outer += [a[i] + [b[j]]]
    return outer

class MultiIndex:

    """A MultiIndex represents a list of indices and holds the
    following data:

        rank    - rank of multiindex
        dims    - a list of dimensions
        indices - list of all possible multiindex values"""

    def __init__(self, dims):
        "Create multi index from given list of dimensions"
        self.rank = len(dims)
        self.dims = listcopy(dims)
        self.indices = build_indices(dims)
        return

    def __repr__(self):
        "Pretty print"
        return "rank = %d dims = %s" % (self.rank, str(self.dims))
