__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-09-29 -- 2007-02-06"
__copyright__ = "Copyright (C) 2004-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

class Index:
    """An Index represents a tensor index. The type of index can be
    either fixed, primary, or secondary as listed below:

      fixed:      tensor rank 0, index is a fixed integer
      primary:    tensor rank 1, index is part of first multiindex (i)
      secondary:  tensor rank 1, index is part of second multiindex (alpha)
      auxiliary:  tensor rank 1, index is part of third multiindex (beta)
      function:   unique index for each Function (coefficient)
      projection: unique index for each projection of a Function (coefficient)
      constant:   unique index for each Constant

      reference tensor auxiliary: auxiliary index local to the reference tensor
      geometry tensor auxiliary:  auxiliary index local to the geometry tensor
      
    The type of index is determined by the arguments to the
    constructor:

      i = Index(0)            creates a given fixed index (0 in this case)
      i = Index("primary")    creates a free primary index
      i = Index("secondary")  creates a free secondary index
      i = Index("function")   creates a Function index
      i = Index("projection") creates a Projection index
      i = Index("constant")   creates a Constant index
      i = Index()             creates a free secondary index (default)

    Note that only primary and secondary indices are used in the
    creation of an element of the algebra. Auxiliary indices are
    detected automatically in a preprocessing of an element of the
    algebra (using the reassign module)."""

    # Avaiable index types
    FIXED       = 0
    PRIMARY     = 1
    SECONDARY   = 2
    AUXILIARY   = 3
    FUNCTION    = 4
    PROJECTION  = 5
    CONSTANT    = 6
    AUXILIARY_0 = 7
    AUXILIARY_G = 8

    def __init__(self, index = "secondary", range = None):
        "Create Index."
        if isinstance(index, Index):
            # Create Index from Index (copy constructor)
            self.index = index.index
            self.type = index.type
            self.range = index.range
        elif isinstance(index, int):
            # Create fixed Index
            self.index = index
            self.type = self.FIXED
            self.range = [index]
        elif index == "primary":
            # Create primary Index
            self.index = next_primary_index()
            self.type = self.PRIMARY
            self.range = range
        elif index == "secondary":
            # Create secondary Index
            self.index = next_secondary_index()
            self.type = self.SECONDARY
            self.range = range
        elif index == "function":
            # Create Function Index
            self.index = next_function_index()
            self.type = self.FUNCTION
            self.range = range
        elif index == "projection":
            # Create Projection Index
            self.index = next_projection_index()
            self.type = self.PROJECTION
            self.range = range
        elif index == "constant":
            # Create Constant Index
            self.index = next_constant_index()
            self.type = self.CONSTANT
            self.range = range # meg: ?
        elif index == "auxiliary":
            # Create auxiliary Index (not possible)
            raise RuntimeError, "Auxiliary indices cannot be created (only modified)."
        elif index == None:
            # Create secondary Index (default)
            self.index = next_secondary_index()
            self.type = self.SECONDARY
            self.range = range
        else:
            raise RuntimeError, "Unknown index type " + str(index)
        return

    def __call__(self, i = [], a = [], b0 = [], b1 = []):
        "Evaluate Index at current index list."
        if self.type == self.FIXED:
            return self.index
        elif self.type == self.PRIMARY:
            if not i:
                raise RuntimeError, "Missing index values for primary indices."
            return i[self.index]
        elif self.type == self.SECONDARY:
            if not a:
                raise RuntimeError, "Missing index values for secondary indices."
            return a[self.index]
        elif self.type == self.AUXILIARY_0:
            if not b0:
                raise RuntimeError, "Missing index values for auxiliary indices."
            return b0[self.index]
        elif self.type == self.AUXILIARY_G:
            if not b1:
                raise RuntimeError, "Missing index values for auxiliary indices."
            return b1[self.index]
        else:
            raise RuntimeError, "Unknown index type " + str(self.type)
        return

    def __cmp__(self, other):
        "Check if Indices are equal."
        # meg: Question, are two indices equal if they have different range?
        if not isinstance(other, Index):
            return -1
        if self.index == other.index and self.type == other.type:
            return 0
        return -1 # Ignore self > other

    def __repr__(self):
        "Print nicely formatted representation of Index."
        offset = "" # meg: Fixme
        if self.range and self.range[0]:
            offset = " + " + str(self.range[0])
        if self.type == self.FIXED:
            return str(self.index)
        elif self.type == self.PRIMARY:
            return "i" + str(self.index) + offset
        elif self.type == self.SECONDARY:
            return "a" + str(self.index) + offset
        elif self.type == self.FUNCTION:
            return str(self.index) + offset
        elif self.type == self.PROJECTION:
            return str(self.index) + offset
        elif self.type == self.CONSTANT:
            return str(self.index) + offset        
        else:
            return "b" + str(self.index) + offset

next_index_0 = 0 # Next available primary index
next_index_1 = 0 # Next available secondary index
next_index_2 = 0 # Next available Function index
next_index_3 = 0 # Next available Projection index
next_index_4 = 0 # Next available Constant index

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

def next_function_index():
    "Return next available Function index."
    global next_index_2
    next_index_2 += 1
    return next_index_2 - 1

def next_projection_index():
    "Return next available Projection index."
    global next_index_3
    next_index_3 += 1
    return next_index_3 - 1

def next_constant_index():
    "Return next available Constant index."
    global next_index_4
    next_index_4 += 1
    return next_index_4 - 1

def reset():
    "Reset all indices."
    global next_index_0, next_index_1, next_index_2, next_index_3, next_index_4
    next_index_0 = 0
    next_index_1 = 0
    next_index_2 = 0
    next_index_3 = 0
    next_index_4 = 0
    return
