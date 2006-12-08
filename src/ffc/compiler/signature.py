"""This module contains utilities for computing signatures of
products. Signatures are used to uniquely identify reference
tensors that may be common to a group of terms."""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2005-09-06 -- 2006-04-04"
__copyright__ = "Copyright (C) 2005-2006 Anders Logg"
__license__  = "GNU GPL Version 2"

# Python modules
from re import sub

# FFC modules
from finiteelement import *
from algebra import *
from tokens import *

def compute_hard_signature(product):
    "Compute hard (unique) signature."
    
    # Create signature for numeric constant
    numeric = "%.15e" % product.numeric
    
    # Create signatures for basis functions
    factors = []
    for v in product.basisfunctions:
        factors += ["{%s;%s;%s;%s,%s}" % \
                    (str(v.element),     \
                     str(v.index),       \
                     str(v.component),   \
                     str(v.derivatives),
                     str(v.restriction))]


    # Sort signatures for basis functions
    factors.sort()

    # Create signature for integral
    integral = str(product.integral)

    # Create signature for product
    return "*".join([numeric] + factors + [integral])

def compute_soft_signature(product):
    "Compute soft (modulo secondary index numbers) signature."

    # Compute hard signature
    signature = compute_hard_signature(product)

    # Ignore secondary index numbers
    return sub('a\d+', 'a', signature)
