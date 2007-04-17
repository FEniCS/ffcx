"""This module implementsq reordering of indices. Reordering of
indices may be necessary in order to factor out common reference
tensors from terms that have the same tensor structure but with
different names of indices."""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2005-09-06 -- 2007-03-05"
__copyright__ = "Copyright (C) 2004-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# Python modules
import sys

# FFC common modules
from ffc.common.exceptions import *
from ffc.common.utils import *
from ffc.common.debug import *

# FFC language modules
from ffc.compiler.language.reassignment import *
from ffc.compiler.language.signature import *

def factorize(monomials):
    """Reorder indices to factorize common reference tensors and
    return a list mapping each term to a previous matching term (if
    any)"""

    debug("Computing factorization...")

    # Create empty factorization
    factorization = [None for i in range(len(monomials))]

    # Compare signatures for pairs of terms
    for i in range(len(monomials) - 1):
        p = monomials[i]
        p_soft = compute_soft_signature(p)
        p_hard = compute_hard_signature(p)

        # Compare term j against term i for i < j
        for j in range(i + 1, len(monomials)):
            
            # Don't factorize against another term if already factorized
            if not factorization[j] == None:
                continue

            # Compute signatures
            q = monomials[j]
            q_soft = compute_soft_signature(q)
            q_hard = compute_hard_signature(q)

            # Check for matching signatures
            if p_hard == q_hard:

                # Group terms if hard signature matches
                factorization[j] = i
                
            elif p_soft == q_soft:

                # Reorder terms if soft signature matches
                monomials[j] = __reorder_indices(p, q, p_hard)
                q = monomials[j]

                # Check that the hard signatures now match
                q_hard = compute_hard_signature(q)
                if not p_hard == q_hard:
                    raise FormError, (form, "Hard signatures don't match after reordering.")

                # Group terms if hard signature matches
                factorization[j] = i

    debug("done")

    return factorization

def __reorder_indices(p, q, p_hard):
    """Reorder secondary indices of monomial q to match the secondary
    indices of monomial p"""

    # Get the number of secondary indices (assuming indices have been
    # previously reassigned to be in the range 0,...,n-1)
    p_max = max_index(p, Index.SECONDARY)
    q_max = max_index(q, Index.SECONDARY)
    if not p_max == q_max:
        raise FormError, ((p, q), "Terms have different index ranges but common soft signature.")
    n = p_max + 1

    # Generate all permutations of indices in the range 0,...,n-1
    for reordering in permutations(range(n)):

        # Copy q and add n to indices (so we can reorder properly)
        q_new = Monomial(q)
        for i in range(n):
            num_changed = reassign_index(q_new, i, i + n, Index.SECONDARY)
            if not num_changed == 1:
                raise FormError, ((p, q), "Not exactly one index modified.")

        # Reorder according to the current reordering
        for i in range(n):
            num_changed = reassign_index(q_new, i + n, reordering[i], Index.SECONDARY)

        # Compare hard signatures for p and q_new
        q_new_hard = compute_hard_signature(q_new)
        if q_new_hard == p_hard:
            return q_new

    raise FormError, ((p, q), "Unable to find a proper reordering of indices.")
