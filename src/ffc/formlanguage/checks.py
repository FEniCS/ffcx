"""This module provides a range of check for a given Sum."""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2006-12-01 -- 2006-12-01"
__copyright__ = "Copyright (C) 2006 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC common modules
from ffc.common.exceptions import *

# FFC compiler modules
from integral import *

def check_sum(sum):
    "Check that given Sum is valid."
    check_integrals(sum)
    check_restrictions(sum)

def check_integrals(sum):
    "Check that all terms have integrals."
    for p in sum.products:
        if not p.integral:
            raise FormError, (p, "Missing integral in term.")

def check_restrictions(sum):
    "Check that all terms are restricted correctly."
    for p in sum.products:
        type = p.integral.type
        for v in p.basisfunctions:
            if type == Integral.CELL:
                if not v.restriction == None:
                    raise FormError, (p, "Integrand may not be restricted in a cell integral.")
            elif type == Integral.EXTERIOR_FACET:
                if not v.restriction == None:
                    raise FormError, (p, "Integrand may not be restricted in an exterior facet integral.")
            elif type == Integral.INTERIOR_FACET:
                if v.restriction == None:
                    raise FormError, (p, "Integrand must be restricted ('+') or ('-') in an interior facet integral.")
