"""This module contains utilities for finding the FiniteElements
defining the test, trial, and function finite element spaces of a
given Sum. All functions assume that all Indices of the given Sum have
already been reassigned."""

__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2005-03-15"
__copyright__ = "Copyright (c) 2005 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC common modules
from ffc.common.debug import *

def find_test(sum):
    "Return the FiniteElement associated with the test function (if any)."
    element = None
    for p in sum.products:
        count = 0
        for v in p.basisfunctions:
            if v.index.type == "primary" and v.index.index == 0:
                count += 1
                if count > 1:
                    raise RuntimeError, "There can only be one test function."
                if element and (not element == v.element):
                    raise RuntimeError, "Test function defined by multiple elements."
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
            if v.index.type == "primary" and v.index.index == 1:
                count += 1
                if count > 1:
                    raise RuntimeError, "There can only be one trial function."
                if element and (not element == v.element):
                    raise RuntimeError, "Trial function defined by multiple elements."
                element = v.element

    if element:
        debug("Finite element of trial space: " + str(element), 0)

    return element

def find_elements(sum, nfunctions):
    "Return a list of FiniteElements associated with all Functions."

    # List of elements used for functions
    elements = [None for j in range(nfunctions)]

    # Iterate over all Coefficients in all Products
    for p in sum.products:
        for c in p.coefficients:
            elements[c.number.index] = c.element

    # Check that we found an element for each function
    for element in elements:
        if not element:
            raise RuntimeError, "Unable to find element for each function."

    if elements:
        debug("Finite elements for functions: " + str(elements), 0)
          
    #return (elements, felement)
    return elements
