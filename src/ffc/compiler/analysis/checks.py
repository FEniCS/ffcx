"""This module provides a number of checks that can be used to check
if a given form is valid."""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2006-12-01 -- 2007-02-06"
__copyright__ = "Copyright (C) 2006-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC common modules
from ffc.common.debug import *
from ffc.common.exceptions import *

# FFC language modules
from ffc.compiler.language.algebra import *
from ffc.compiler.language.integral import *
from ffc.compiler.language.indexcall import *

def check_form(form):
    "Check that the form is valid"
    debug("Checking validity of form...")
    check_type(form)
    check_integrals(form)
    check_restrictions(form)
    check_completeness(form)
    debug("ok")

def check_type(form):
    "Check that the form is a Form"
    if not isinstance(form, Form):
        raise FormError, (p, "Not a form.")

def check_integrals(form):
    "Check that all terms have integrals"
    for p in form.monomials:
        if not p.integral:
            raise FormError, (p, "Missing integral in term.")

def check_restrictions(form):
    "Check that all terms are restricted correctly"
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

def check_completeness(form):
    "Check that each secondary index appears exactly twice in each term"
    for m in form.monomials:
        aindices = []
        bindices = []
        index_call(m, index_add, [aindices, Index.SECONDARY])
        index_call(m, index_add, [bindices, Index.AUXILIARY])
        if not __check_completeness(aindices) or not __check_completeness(bindices):
            raise FormError, (m, "Index does not appear exactly twice in term.")

def __check_completeness(indices):
    "Check that each index in the list appear exactly twice"
    for i in range(len(indices)):
        count = 0
        for j in range(len(indices)):
            if indices[i] == indices[j]:
                count += 1
        if not count == 2:
            return False
    return True
