"""This module contains utilities for finding the FiniteElements
defining the test, trial, and function finite element spaces of a
given Sum. All functions assume that all Indices of the given Sum have
already been reassigned."""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2005-03-15 -- 2007-01-23"
__copyright__ = "Copyright (C) 2005-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# Modified by Garth N. Wells 2006

# Python modules
import sys
import numpy

# FFC common modules
sys.path.append("../../")
from ffc.common.debug import *
from ffc.common.exceptions import *
from ffc.common.constants import *

# FFC form language modules
from ffc.formlang.index import *

# FFC compiler modules
from declaration import *

def find_test(sum):
    "Return the FiniteElement associated with the test function (if any)."
    element = None
    for p in sum.products:
        count = 0
        for v in p.basisfunctions:
            if v.index.type == Index.PRIMARY and v.index.index == 0:
                count += 1
                if count > 1:
                    raise FormError, (sum , "There can only be one test function.")
                if element and (not element == v.element):
                    raise FormError, (sum, "Test function defined by multiple elements.")
                element = v.element

    if element:
        debug("Finite element of test space:  " + str(element), 0)

    return element

def find_trial(sum):
    "Return the FiniteElement associated with the trial function."
    element = None
    for p in sum.products:
        count = 0
        for v in p.basisfunctions:
            if v.index.type == Index.PRIMARY and v.index.index == 1:
                count += 1
                if count > 1:
                    raise FormError, (sum, "There can only be one trial function.")
                if element and (not element == v.element):
                    raise FormError, (sum, "Trial function defined by multiple elements.")
                element = v.element

    if element:
        debug("Finite element of trial space: " + str(element), 0)

    return element

def find_elements(sum, nfunctions):
    """Return a list of FiniteElements associated with the (original)
    function spaces of the Functions appearing in the form."""

    # List of elements used for functions
    elements = [None for j in range(nfunctions)]

    # Iterate over all Coefficients in all Products
    for p in sum.products:
        for c in p.coefficients:
            elements[c.n0.index] = c.e0

    # Check that we found an element for each function
    for element in elements:
        if not element:
            raise FormError, (sum, "Unable to find element for each function.")

    if elements:
        debug("Finite elements for functions: " + str(elements), 0)
          
    return elements

def find_projections(sum, nprojections):
    """Return a list of tuples (n0, n1, e0, e1, P) defining the
    projections of all Functions appearing in the form."""

    # List of projections used for functions
    projections = [None for j in range(nprojections)]

    # Iterate over all Coefficients in all Products
    for p in sum.products:
        for c in p.coefficients:
            projections[c.n1.index] = (c.n0.index, c.n1.index, c.e0, c.e1, c.P)

    # Check that we found an element for each projection
    for projection in projections:
        if not projection:
            raise FormError, (sum, "Unable to find a projection for each function.")

    #if projections:
    #    debug("Projections for functions: " + str(projections), 0)
          
    return projections
