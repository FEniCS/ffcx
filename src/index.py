__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-09-29"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# Python modules
from Numeric import *

next_index_0 = 0
next_index_1 = 0

def next_primary_index():
    "Return next available primary index."
    global next_index_0
    next_index_0 += 1
    return next_index_0 - 1

def next_secondary_index():
    "Return next available secondary index."
    global next_index_1
    next_index_1 += 1
    return next_index_1 - 1

def reset():
    "Reset all indices."
    global next_index_0, next_index_1
    next_index_0 = 0
    next_index_1 = 0
    return

def build_indices(dims):
    """Create a list of all index combinations matching the given list
    of index dimensions.  Someone please tell me if there is a better
    way to iterate over a multi-dimensional array. The list of indices
    is constucted by iterating over all integer numbers that can be
    represented with the base of each position determined by the given
    dimension for that index."""
    if not dims:
        return []
    current = zeros(len(dims))
    indices = []
    posvalue = [1] + list(cumproduct(dims)[:-1])
    for i in range(product(dims)):
        sum = i
        for pos in range(len(dims)):
            current[pos] = (i / posvalue[pos]) % dims[pos]
            sum -= current[pos] * posvalue[pos]
        indices += [list(current)]
    return indices

class Index:
    
    """An Index represents a tensor index. The type of index can be
    either fixed, primary, or secondary as listed below:

      fixed:     tensor rank 0, index is a fixed integer
      primary:   tensor rank 1, index is part of first multiindex (i)
      secondary: tensor rank 1, index is part of second multiindex (alpha)
      auxiliary: tensor rank 1, index does not appear inside the integral
      
    The type of index is determined by the arguments to the
    constructor:

      i = Index(0)           creates a given fixed index (0 in this case)
      i = Index("primary")   creates a free primary index
      i = Index("secondary") creates a free secondary index
      i = Index("auxiliary") creates a free auxiliary index
      i = Index()            creates a free secondary index

    Note that only primary and secondary indices are used in the
    creation of an element of the algebra. Auxiliary indices are
    detected automatically in a preprocessing of an element of the
    algebra when creating a Form (using the reassign module)."""

    def __init__(self, index = "secondary"):
        "Create Index."
        if isinstance(index, Index):
            self.index = index.index
            self.type = index.type
        elif isinstance(index, int):
            self.index = index
            self.type = "fixed"
        elif index == "primary":
            self.index = next_primary_index()
            self.type = "primary"
        elif index == "secondary":
            self.index = next_secondary_index()
            self.type = "secondary"
        elif index == "auxiliary":
            raise RuntimeError, "Auxiliary indices cannot be created (only modified)."
        elif index == None:
            self.index = next_secondary_index()
            self.type = "secondary"
        else:
            raise RuntimeError, "Unknown index type " + str(index)
        return

    def __call__(self, indices0, indices1, indices2):
        "Evaluate Index at current index list."
        if self.type == "primary":
            if not indices0:
                raise RuntimeError, "Missing index values for primary indices."
            return indices0[self.index]
        elif self.type == "secondary":
            if not indices1:
                raise RuntimeError, "Missing index values for secondary indices."
            return indices1[self.index]
        elif self.type == "auxiliary":
            if not indices2:
                raise RuntimeError, "Missing index values for auxiliary indices."
            return indices2[self.index]
        else:
            raise RuntimeError, "Uknown index type " + str(self.type)
        return
        
    def __repr__(self):
        "Print nicely formatted representation of Index."
        if self.type == "fixed":
            return str(self.index)
        elif self.type == "primary":
            return "i" + str(self.index)
        elif self.type == "secondary":
            return "a" + str(self.index)
        elif self.type == "auxiliary":
            return "b" + str(self.index)
        else:
            raise RuntimeError, "Unknown index type " + str(index)
        return
