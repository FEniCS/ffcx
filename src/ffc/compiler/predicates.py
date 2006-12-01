"""This module provides predicates for forms."""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2006-12-01 -- 2006-12-01"
__copyright__ = "Copyright (C) 2006 Anders Logg"
__license__  = "GNU GPL Version 2"

def has_integrals(sum):
    "Check if all terms have integrals."
    for p in sum.products:
        if not p.integral:
            return False
    return True
