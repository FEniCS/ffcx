"""This module contains utilities for finding the FiniteElements
defining the test, trial, and function finite element spaces of a
given Sum. All functions assume that all Indices of the given Sum have
already been reassigned."""

__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2005-03-15 -- 2005-11-08"
__copyright__ = "Copyright (c) 2005 Anders Logg"
__license__  = "GNU GPL Version 2"

# Python modules
import sys
import Numeric

# FFC common modules
sys.path.append("../../")
from ffc.common.debug import *
from ffc.common.exceptions import *
from ffc.common.constants import *

# FFC compiler modules
from declaration import *

def find_test(sum):
    "Return the FiniteElement associated with the test function (if any)."
    element = None
    for p in sum.products:
        count = 0
        for v in p.basisfunctions:
            if v.index.type == "primary" and v.index.index == 0:
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
            if v.index.type == "primary" and v.index.index == 1:
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
    """Return a tuple (elements, projections) consisting of a list of
    FiniteElements associated with all Functions and a list of
    corresponding projections, where projections[i] is the projection
    matrix P from elements[i] to the space used in the representation
    of the form."""

    # List of elements used for functions
    elements = [None for j in range(nfunctions)]

    # List of corresponding projections
    projections = [None for j in range(nfunctions)]

    # Iterate over all Coefficients in all Products
    for p in sum.products:
        for c in p.coefficients:
            elements[c.number.index] = c.e0 or c.element
            projections[c.number.index] = c.projection

    # Check that we found an element for each function
    for element in elements:
        if not element:
            raise FormError, (sum, "Unable to find element for each function.")

    if elements:
        debug("Finite elements for functions: " + str(elements), 0)
          
    return (elements, projections)

def compute_coefficients(elements, projections, format):
    "Precompute declarations of coefficients according to given format."

    declarations = []

    for j in range(len(elements)):

        element = elements[j]
        P = projections[j]

        if P == None:
            # No projection, just copy the values
            for k in range(element.spacedim()):
                name = format.format["coefficient"](j, k)
                value = format.format["coefficient table"](j, k)
                declarations += [Declaration(name, value)]
        else:
            # Compute projection
            (m, n) = Numeric.shape(P)
            for k in range(m):
                terms = []
                for l in range(n):
                    if abs(P[k][l] < FFC_EPSILON):
                        continue
                    cl = format.format["coefficient table"](j, l)
                    if abs(P[k][l] - 1.0) < FFC_EPSILON:
                        terms += [cl]
                    else:
                        Pkl = format.format["floating point"](P[k][l])
                        terms += [format.format["multiplication"]([Pkl, cl])]
                name = format.format["coefficient"](j, k)
                value = format.format["sum"](terms)
                declarations += [Declaration(name, value)]

    return declarations
