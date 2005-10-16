"""This module contains utilities for reordering of indices.
Reordering of indices may be necessary in order to factor out common
reference tensors from terms that have the same tensor structure but
with different names of indices."""

__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2005-09-06 -- 2005-09-14"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# Python modules
import sys

# FFC common modules
sys.path.append("../../")
from ffc.common.debug import *
from ffc.common.exceptions import *
from ffc.common.util import permutations

# FFC compiler modules
from reassign import *
from signature import *

def reorder_indices(sum):
    """Reorder indices to find terms with common reference tensors and
    compute factorization, mapping terms with already computed
    reference tensors to the already computed matching term."""

    if not isinstance(sum, Sum):
        raise FormError, (sum, "Indices can only be reordered for Sums.")

    # Create empty factorization
    factorization = [None for i in range(len(sum.products))]

    # Compare signatures for pairs of terms
    for i in range(len(sum.products) - 1):

        p = sum.products[i]
        p_soft = compute_soft_signature(p)
        p_hard = compute_hard_signature(p)
        debug("Soft signature: " + p_soft, 1)
        debug("Hard signature: " + p_hard, 1)

        # Compare term i against term j for j > i
        for j in range(i + 1, len(sum.products)):

            # Don't factorize against another term if already factorized
            if not factorization[j] == None:
                continue

            # Compute signatures
            q = sum.products[j]
            q_soft = compute_soft_signature(q)
            q_hard = compute_hard_signature(q)

            if p_hard == q_hard:
                print "Hard signatures match for terms %d and %d, factorizing" % (i, j)
                # Group terms if hard signature matches
                factorization[j] = i
            elif p_soft == q_soft:
                print "Soft signatures match for terms %d and %d, reordering and factorizing" % (i, j)
                # Reorder terms if soft signature matches
                sum.products[j] = __reorder_indices(p, q, p_hard)
                q = sum.products[j]
                # Check that the hard signatures now match
                q_hard = compute_hard_signature(q)
                if not p_hard == q_hard:
                    raise FormError, (sum, "Hard signatures don't match after reordering.")
                # Group terms if hard signature matches
                factorization[j] = i

    return factorization

def __reorder_indices(p, q, p_hard):
    """Reorder secondary indices of Product q to match the secondary
    indices of Product p."""

    # Get the number of secondary indices (assuming indices have been
    # previously reassigned to be in the range 0,...,n-1)
    p_max = max_index(p, "secondary")
    q_max = max_index(q, "secondary")
    if not p_max == q_max:
        raise FormError, ((p, q), "Terms have different index ranges but common soft signature.")
    n = p_max + 1

    # Generate all permutations of indices in the range 0,...,n-1
    for reordering in permutations(range(n)):
        #print "Reordering indices according to " + str(reordering)
        # Copy q and add n to indices (so we can reorder properly)
        q_new = Product(q)
        for i in range(n):
            num_changed = reassign_index(q_new, i, i + n, "secondary")
            if not num_changed == 1:
                raise FormError, ((p, q), "Not exactly one index modified.")
        # Reorder according to the current reordering
        for i in range(n):
            num_changed = reassign_index(q_new, i + n, reordering[i], "secondary")
        # Compare hard signatures for p and q_new
        q_new_hard = compute_hard_signature(q_new)
        if q_new_hard == p_hard:
            return q_new

    raise Form, ((p, q), "Unable to find a proper reordering of indices.")
