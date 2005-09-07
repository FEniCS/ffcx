"""This module contains utilities for reordering of indices.
Reordering of indices may be necessary in order to factor out common
reference tensors from terms that have the same tensor structure but
with different names of indices."""

__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2005-09-06"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC common modules
from ffc.common.util import permutations

# FFC compiler modules
from reassign import *
from signature import *

def reorder_indices(sum):
    "Reorder indices to find terms with common reference tensors."
    if not isinstance(sum, Sum):
        raise RuntimeError, "Indices can only be reordered for Sums."

    # Compare signatures for pairs of terms
    for i in range(len(sum.products) - 1):

        p = sum.products[i]
        p_soft = compute_soft_signature(p)
        p_hard = compute_hard_signature(p)

        # Compare term i against term j for j > i
        for j in range(i + 1, len(sum.products)):

            q = sum.products[j]
            q_soft = compute_soft_signature(q)
            q_hard = compute_hard_signature(q)

            # Group terms if hard signature matches
            if p_hard == q_hard:
                print "Hard signatures match for terms %d and %d" % (i, j)

            # Reorder terms if soft signature matches, then group
            if p_soft == q_soft:
                __reorder_indices(p, q)
                print "Soft signatures match for terms %d and %d" % (i, j)

    return

def __reorder_indices(p, q):
    """Reorder secondary indices of Product q to match the secondary
    indices of Product p."""

    # Get the number of secondary indices (assuming indices have been
    # previously reassigned to be in the range 0,...,n-1)
    p_max = max_index(p, "secondary")
    q_max = max_index(q, "secondary")
    if not p_max == q_max:
        raise RuntimeError, \
              "Terms have different index ranges but common soft signature"
    n = p_max + 1

    # Generate all permutations of indices in the range 0,...,n-1
    for reordering in permutations(range(n)):
        print "Reordering indices according to " + str(reordering)
