"""This module contains utilities for reordering of indices.
Reordering of indices may be necessary in order to factor out common
reference tensors from terms that have the same tensor structure but
with different names of indices."""

__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2005-09-06"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC modules
from signature import *

def reorder_indices(sum):
    "Reorder indices to find terms with common reference tensors."
    if not isinstance(sum, Sum):
        raise RuntimeError, "Indices can only be reordered for Sums."

    for p in sum.products:
        soft = compute_soft_signature(p)
        hard = compute_hard_signature(p)

        print "Soft signature: " + soft
        print "Hard signature: " + hard

    return
