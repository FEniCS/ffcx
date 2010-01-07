__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2005-02-04"
__copyright__ = "Copyright (C) 2005-2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian B. Oelgaard, 2009
# Last changed: 2010-01-07

# Python modules.
import operator

# FFC modules.
from log import error

def product(sequence):
    "Return the product of all elements in a sequence."
    # Copied from UFL
    return reduce(operator.__mul__, sequence, 1)

def pick_first(values):
    "Check that all values are equal and return the value."
    if not values[:-1] == values[1:]:
        error("Values differ: " + str(values))
    return values[0]

def listcopy(l):
    """Create a copy of the list, calling the copy constructor on each
    object in the list (problems when using copy.deepcopy)."""
    if not l:
        return []
    else:
        return [object.__class__(object) for object in l]

def compute_permutations(k, n, skip = []):
   """Compute all permutations of k elements from (0, n) in rising order.
   Any elements that are contained in the list skip are not included."""
   if k == 1:
       return [(i,) for i in range(n) if not i in skip]
   pp = compute_permutations(k - 1, n, skip)
   permutations = []
   for i in range(n):
       if i in skip:
           continue
       for p in pp:
           if i < p[0]:
               permutations += [(i, ) + p]
   return permutations
