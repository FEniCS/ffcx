"""This module provides a range of check for a given Form."""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2006-12-01 -- 2007-01-23"
__copyright__ = "Copyright (C) 2006-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC common modules
from ffc.common.exceptions import *

# FFC compiler modules
from integral import *

def check_form(form):
    "Check that given Form is valid."
    check_integrals(form)
    check_restrictions(form)

def check_integrals(form):
    "Check that all terms have integrals."
    for p in form.monomials:
        if not p.integral:
            raise FormError, (p, "Missing integral in term.")

def check_restrictions(form):
    "Check that all terms are restricted correctly."
    for p in form.monomials:
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
