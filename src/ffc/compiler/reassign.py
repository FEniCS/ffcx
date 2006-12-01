"""This module contains utilities for reassignment of indices.
Reassignment is necessary since when an element of the algebra is
created, new indices are created with unique numbers not necessarily
starting at 0. This makes it simpler to implement the algebra, since
the algebra doesn't need to keep track of indices."""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-10-13 -- 2006-12-01"
__copyright__ = "Copyright (C) 2004-2006 Anders Logg"
__license__  = "GNU GPL Version 2"

# Python modules
import sys

# FFC common modules
sys.path.append("../../")
from ffc.common.exceptions import *

# FFC compiler modules
import algebra
from index import *

def reassign_indices(sum):
    """Reassign indices of Sum with all indices starting at 0. Modify
    secondary indices not appearing in both the reference and geometry
    tensors to auxiliary indices."""

    # Check that we got a Sum
    if not isinstance(sum, algebra.Sum):
        raise RuntimeError, "Can only reassign indices for Sum."

    # Check for completeness
    if not iscomplete(sum):
        raise FormError, (sum, "Given form is not complete.")

    # Modify secondary indices to auxiliary indices
    [__create_auxiliary(p) for p in sum.products]

    # Reassign primary indices (global for the Sum)
    __reassign(sum, Index.PRIMARY)

    # Reassign Function indices (global for the Sum)
    __reassign(sum, Index.FUNCTION)

    # Reassign Projection indices (global for the Sum)
    __reassign(sum, Index.PROJECTION)

    # Reassign Constant indices (global for the Sum)
    __reassign(sum, Index.CONSTANT)

    # Reassign secondary indices (global for each Product)
    [__reassign(p, Index.SECONDARY) for p in sum.products]

    # Reassign reference tensor auxiliary indices (global for each Product)
    [__reassign(p, Index.AUXILIARY_0) for p in sum.products]

    # Reassign geometry tensor auxiliary indices (global for each Product)
    [__reassign(p, Index.AUXILIARY_G) for p in sum.products]

    # Check for completeness again
    if not iscomplete(sum):
        raise FormError, (sum, "Given form is not complete.")
    
    return

def iscomplete(object):
    """Check if given Product or Sum is complete with respect to
    secondary and auxiliary Indices, that is, each such Index appear
    exactly twice in each Product."""
    if isinstance(object, algebra.Sum):
        # Check that each Product is complete
        for p in object.products:
            if not iscomplete(p):
                return False
        return True
    elif isinstance (object, algebra.Product):
        # Get secondary and auxiliary Indices
        aindices = []
        bindices = []
        object.indexcall(__index_add, [aindices, Index.SECONDARY])
        object.indexcall(__index_add, [bindices, Index.AUXILIARY])
        return __check_completeness(aindices) and __check_completeness(bindices)
    else:
        raise RuntimeError, "Unsupported type."

def reassign_complete(product, type):
    """Reassign complete secondary and auxiliary Index pairs of given
    type for the given Product, so that they don't collide with
    Indices in other Products (that may get multiplied with current
    Product)."""
    # Get indices
    indices = []
    product.indexcall(__index_add, [indices, type])
    # Find complete pairs
    pairs = []
    for i in range(len(indices)):
        count = 0
        matching = None
        for j in range(i + 1, len(indices)):
            if indices[i] == indices[j]:
                count += 1
                matching = indices[j]
        if count == 1:
            pairs += [(indices[i], matching)]
    # Reassign complete pairs
    for pair in pairs:
        # Get next available index
        if type == Index.SECONDARY:
            new_index = next_secondary_index()
        elif type == Index.AUXILIARY:
            new_index = new_auxiliary_index()
        else:
            raise RuntimeError, "Reassigning for wrong Index type."
        # Reassign for current pair
        pair[0].index = new_index
        pair[1].index = new_index

def reassign_index(object, iold, inew, type):
    """Change value of index from iold to inew for given object,
    and return the number of indices changed."""
    increment = [0]
    object.indexcall(__index_reassign, (iold, inew, type, increment))
    return increment[0]
    
def min_index(object, type):
    "Compute minimum index of given type for given object."
    indices = []
    object.indexcall(__index_add_value, [indices, type])
    return min([sys.maxint] + indices)

def max_index(object, type):
    "Compute maximum index of given type for given object."
    indices = []
    object.indexcall(__index_add_value, [indices, type])
    return max([-1] + indices)

def __reassign(object, type):
    "Reassign all indices of given type for object from iold to inew."
    imin = min_index(object, type)
    imax = max_index(object, type)
    inew = 0 # Increase when index is found
    for iold in range(imin, imax + 1):
        increment = [0]
        object.indexcall(__index_reassign, [iold, inew, type, increment])
        inew += increment[0]
    imin = min_index(object, type)
    imax = max_index(object, type)
    if 0 < imin < sys.maxint:
        raise RuntimeError, "Failed to reassign indices."

def __have_index(object, index):
    "Check if object contains given index."
    indices = []
    object.indexcall(__index_add_value, [indices, index.type])
    return index.index in indices

def __create_auxiliary(product):
    """Modify secondary indices not appearing in both the reference
    tensor and geometry tensor to auxiliary indices."""

    # Build list of reference tensor secondary indices
    ir = []
    [v.indexcall(__index_add_value, [ir, Index.SECONDARY]) for v in product.basisfunctions]
    # Build list of geometry tensor secondary indices
    ig = []
    [c.indexcall(__index_add_value, [ig, Index.SECONDARY]) for c in product.coefficients]
    [t.indexcall(__index_add_value, [ig, Index.SECONDARY]) for t in product.transforms]
    # Build lists of indices that should be modified
    irm = [i for i in ir if i not in ig]
    igm = [i for i in ig if i not in ir]
    # Modify indices of reference tensor
    for v in product.basisfunctions:
        v.indexcall(__index_modify, (irm, Index.AUXILIARY_0))
    # Modify indices of geometry tensor
    for c in product.coefficients:
        c.indexcall(__index_modify, (igm, Index.AUXILIARY_G))
    for t in product.transforms:
        t.indexcall(__index_modify, (igm, Index.AUXILIARY_G))
    return

def __index_add(index, args):
    "Add Index to list if Index is of given type."
    indices = args[0]
    type = args[1]
    if index.type == type:
        indices += [index]
    return

def __index_add_value(index, args):
    "Add Index number to list if index is of given type."
    indices = args[0]
    type = args[1]
    if index.type == type:
        indices += [index.index]
    return

def __index_modify(index, args):
    "Modify secondary index to index of given type if it appears in the list."
    indices = args[0]
    type = args[1]
    if index.index in indices and index.type == Index.SECONDARY:
        index.type = type
    return

def __index_reassign(index, args):
    "Reassign index from old value to new value."
    iold = args[0]
    inew = args[1]
    type = args[2]
    increment = args[3]
    if index.index == iold and index.type == type:
        index.index = inew
        increment[0] = 1
    return

def __check_completeness(indices):
    "Check that each Index in the list appear exactly twice."
    for i in range(len(indices)):
        count = 0
        for j in range(len(indices)):
            if indices[i] == indices[j]:
                count += 1
        if not count == 2:
            print "Index %s found in %d position(s), but must appear in exactly two positions." % \
                  (str(indices[i]), count)
            return False
    return True

if __name__ == "__main__":

    print "Testing reassignment"
    print "--------------------"
    
    element = FiniteElement("Lagrange", "triangle", 1)
    
    u = BasisFunction(element)
    v = BasisFunction(element)
    i = Index()

    w1 = u.dx(i)*v.dx(i) + u*v
    w2 = u.dx(i)*v.dx(i) + u*v

    print w1
    print w2
    
    reassign_indices(w1)
    reassign_indices(w2)

    print w1
    print w2
