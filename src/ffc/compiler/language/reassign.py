"""This module provides reassignment of indices for forms in such a
way that each unique index is given a number from 0 to n - 1 where n
is the number of indices."""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-10-13 -- 2007-02-06"
__copyright__ = "Copyright (C) 2004-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# Python modules
import sys

# FFC common modules
from ffc.common.debug import *
from ffc.common.exceptions import *

# FFC language modules
from index import *
from indexcall import *

def reassign_indices(form):
    """Reassign indices of Form with all indices starting at 0 and
    modify secondary indices not appearing in both the reference and
    geometry tensors to auxiliary indices"""

    debug("Reassigning form indices...")

    # Modify secondary indices to auxiliary indices
    [__create_auxiliary(m) for m in form.monomials]

    print form

    # Reassign primary indices (global for the Form)
    __reassign(form, Index.PRIMARY)

    # Reassign Function indices (global for the Form)
    __reassign(form, Index.FUNCTION)

    # Reassign Projection indices (global for the Form)
    __reassign(form, Index.PROJECTION)

    # Reassign Constant indices (global for the Form)
    __reassign(form, Index.CONSTANT)

    # Reassign secondary indices (global for each Monomial)
    [__reassign(p, Index.SECONDARY) for p in form.monomials]

    # Reassign reference tensor auxiliary indices (global for each Monomial)
    [__reassign(p, Index.AUXILIARY_0) for p in form.monomials]

    # Reassign geometry tensor auxiliary indices (global for each Monomial)
    [__reassign(p, Index.AUXILIARY_G) for p in form.monomials]

    debug("done")

def reassign_complete(monomial, type):
    """Reassign complete secondary and auxiliary index pairs of given
    type for the given monomial, so that they don't collide with
    indices in other monomials (that may get multiplied with current
    monomial)"""
    
    # Get indices
    indices = []
    index_call(monomial, index_add, [indices, type])
    
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
    and return the number of indices changed"""
    increment = [0]
    index_call(object, index_reassign, (iold, inew, type, increment))
    return increment[0]
    
def min_index(object, type):
    "Compute minimum index of given type for given object"
    indices = []
    index_call(object, index_add_value, [indices, type])
    return min([sys.maxint] + indices)

def max_index(object, type):
    "Compute maximum index of given type for given object"
    indices = []
    index_call(object, index_add_value, [indices, type])
    return max([-1] + indices)

def __reassign(object, type):
    "Reassign all indices of given type for object from iold to inew"
    imin = min_index(object, type)
    imax = max_index(object, type)
    inew = 0 # Increase when index is found
    for iold in range(imin, imax + 1):
        increment = [0]
        index_call(object, index_reassign, [iold, inew, type, increment])
        inew += increment[0]
    imin = min_index(object, type)
    imax = max_index(object, type)
    if 0 < imin < sys.maxint:
        raise RuntimeError, "Failed to reassign indices."

def __have_index(object, index):
    "Check if object contains given index"
    indices = []
    index_call(object, index_add_value, [indices, index.type])
    return index.index in indices

def __create_auxiliary(monomial):
    """Modify secondary indices not appearing in both the reference
    tensor and geometry tensor to auxiliary indices"""

    # Build list of reference tensor secondary indices
    ir = []
    [index_call(v, index_add_value, [ir, Index.SECONDARY]) for v in monomial.basisfunctions]
    
    # Build list of geometry tensor secondary indices
    ig = []
    [index_call(c, index_add_value, [ig, Index.SECONDARY]) for c in monomial.coefficients]
    [index_call(t, index_add_value, [ig, Index.SECONDARY]) for t in monomial.transforms]
    
    # Build lists of indices that should be modified
    irm = [i for i in ir if i not in ig]
    igm = [i for i in ig if i not in ir]
    
    # Modify indices of reference tensor
    for v in monomial.basisfunctions:
        index_call(v, index_modify, (irm, Index.SECONDARY, Index.AUXILIARY_0))

    # Modify indices of geometry tensor
    for c in monomial.coefficients:
        index_call(c, index_modify, (igm, Index.SECONDARY, Index.AUXILIARY_G))
    for t in monomial.transforms:
        index_call(t, index_modify, (igm, Index.SECONDARY, Index.AUXILIARY_G))
    return
