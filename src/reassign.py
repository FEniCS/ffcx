"""This module contains utilities for reassignment of indices and
computation of ranks. Reassignment is necessary since when an Element
of the algebra is created, new indices are created with unique numbers
not necessarily starting at 0. This makes it simpler to implement the
algebra, since the algebra doesn't need to keep track of ranks and
reassign indices."""

__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-10-13"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# Python modules
from sys import maxint

# FFC modules
from algebra import *

def reassign_indices(sum):
    """Reassign indices of Sum with all indices starting at 0. Modify
    secondary indices to auxiliary indices if they don't appear inside
    the integral."""

    # Check that we got a Sum
    if not isinstance(sum, Sum):
        raise RuntimeError, "Can only reassign indices for Sum."
    s = Sum(sum)


    # Modify secondary indices to auxiliary indices
    for i in range(len(s.products)):
        s.products[i] = __create_auxiliary(s.products[i])

    # Reassign primary indices (global for the Sum)
    i0min = min_sum(s, "primary")
    i0max = max_sum(s, "primary")
    i0 = 0
    for i in range(i0min, i0max + 1):
        (s, increment) = __reassign_sum(s, i, i0, "primary")
        i0 += increment # Increases when the index is found
    
    # Reassign secondary indices (global for each Product)
    for j in range(len(s.products)):
        p = s.products[j]
        i1min = min_product(p, "secondary")
        i1max = max_product(p, "secondary")
        i1 = 0
        for i in range(i1min, i1max + 1):
            (p, increment) = __reassign_product(p, i, i1, "secondary")
            i1 += increment #Increases when the index is found
        s.products[j] = p

    # Reassign secondary indices (global for each Product)
    for j in range(len(s.products)):
        p = s.products[j]
        i1min = min_product(p, "auxiliary")
        i1max = max_product(p, "auxiliary")
        i1 = 0
        for i in range(i1min, i1max + 1):
            (p, increment) = __reassign_product(p, i, i1, "auxiliary")
            i1 += increment #Increases when the index is found
        s.products[j] = p

    # Check that indices start at 0 (can be -1 if missing)
    if min_sum(s, "primary") > 0 or min_sum(s, "secondary") or min_sum(s, "auxiliary"):
        raise RuntimeError, "Failed to reassign indices."
    else:
        print "Reassigned indices ok."

    return s

def __reassign_sum(sum, iold, inew, type):
    "Reassign indices of Sum from iold to inew."
    s = Sum(sum)
    result = [__reassign_product(p, iold, inew, type) for p in s.products]
    s.products = [r[0] for r in result]
    increment = max([r[1] for r in result])
    return (s, increment)

def __reassign_product(product, iold, inew, type):
    "Reassign indices of Product from iold to inew."
    increment = 0 # Change to 1 if we find index
    p = Product(product)
    
    # Reassign indices of Factors
    result = [__reassign_factor(f, iold, inew, type) for f in p.factors]
    p.factors = [r[0] for r in result]
    increment = max([increment] + [r[1] for r in result])

    # Reassign indices of Transforms
    result = [__reassign_transform(t, iold, inew, type) for t in p.transforms]
    p.transforms = [r[0] for r in result]
    increment = max([increment] + [r[1] for r in result])

    # Reassign indices of Coefficients
    result = [__reassign_coefficient(c, iold, inew, type) for c in p.coefficients]
    p.coefficients = [r[0] for r in result]
    increment = max([increment] + [r[1] for r in result])

    return (p, increment)

def __reassign_factor(factor, iold, inew, type):
    "Reassign indices of Factor from iold to inew."
    increment = 0 # Change to 1 if we find index
    f = Factor(factor)

    # Reassign index of BasisFunction
    i = f.basisfunction.index
    if i.type == "type" and i.index == iold:
        f.basisfunction.index.index = inew
        increment = 1

    # Reassign indices of Derivatives
    result = [__reassign_derivative(d, iold, inew, type) for d in f.derivatives]
    f.derivatives = [r[0] for r in result]
    increment = max([increment] + [r[1] for r in result])

    # Check that the list of Transforms is empty
    if f.transforms:
        raise RuntimeError, "Non-empty list of Transforms for Factor."

    return (f, increment)

def __reassign_transform(transform, iold, inew, type):
    "Reassign indices of Transform from iold to inew."
    increment = 0 # Change to 1 if we find index
    t = Transform(transform)

    # Check first index
    if t.index0.type == type and t.index0.index == iold:
        t.index0.index = inew
        increment = 1

    # Check second index
    if t.index1.type == type and t.index1.index == iold:
        t.index1.index = inew
        increment = 1

    return (t, increment)

def __reassign_derivative(derivative, iold, inew, type):
    "Reassign index of Derivative from iold to inew."
    increment = 0 # Change to 1 if we find index
    d = Derivative(derivative)
    if d.index.type == type and d.index.index == iold:
        d.index.index = inew
        increment = 1
    return (d, increment)

def __reassign_coefficient(coefficient, iold, inew, type):
    "Reassign index of Coefficient from iold to inew."
    increment = 0 # Change to 1 if we find index
    # FIXME: Do something for the Coefficient
    return (coefficient, increment)

def min_sum(sum, type):
    "Compute minimum index of Sum."
    return min([min_product(p, type) for p in sum.products])

def max_sum(sum, type):
    "Compute maximum index of Sum."
    return max([max_product(p, type) for p in sum.products])

def min_product(product, type):
    "Compute minimum index of Product."
    imin = min([__min_factor(f, type) for f in product.factors])
    imin = min([imin] + [t.index0.index for t in product.transforms if t.index0.type == type])
    imin = min([imin] + [t.index1.index for t in product.transforms if t.index1.type == type])
    # FIXME: imin = min([imin] + [c something for c in product.coefficients])
    return imin

def max_product(product, type):
    "Compute maximum index of Product."
    imax = max([__max_factor(f, type) for f in product.factors])
    imax = max([imax] + [t.index0.index for t in product.transforms if t.index0.type == type])
    imax = max([imax] + [t.index1.index for t in product.transforms if t.index1.type == type])
    # FIXME: imax = max([imax] + [c something for c in product.coefficients])
    return imax

def __min_factor(factor, type):
    "Compute minimum index of Factor."
    imin = maxint
    if factor.basisfunction.index.type == type:
        imin = min(imin, factor.basisfunction.index.index)
    imin = min([imin] + [d.index.index for d in factor.derivatives if d.index.type == type])
    # Check that the list of Transforms is empty
    if factor.transforms:
        raise RuntimeError, "Non-empty list of Transforms for Factor."
    return imin

def __max_factor(factor, type):
    "Compute maximum index of Factor."
    imax = -1
    if factor.basisfunction.index.type == type:
        imax = max(imax, factor.basisfunction.index.index)
    imax = max([imax] + [d.index.index for d in factor.derivatives if d.index.type == type])
    # Check that the list of Transforms is empty
    if factor.transforms:
        raise RuntimeError, "Non-empty list of Transforms for Factor."
    return imax

def dims_product(product, r, type):
    """Compute dimensions for the tensor represented by the
    Product. This method involves some searching, but it
    shouldn't take too much time."""
    dims = []
    for i in range(r):
        (dim, found) = dim_product(product, i, type)
        if found:
            dims = dims + [dim]
        else:
            raise RuntimeError, "Unable to find primary index " + str(i)
    return dims

def dim_product(product, i, type):
    """Try to find dimension number i of the given type. If found, the
    dimension is returned as (dim, True). Otherwise (0, False) is
    returned."""

    # Check Factors
    for f in product.factors:
        (dim, found) = dim_factor(f, i, type)
        if found:
            return (dim, True)

    # Check Transforms, assuming that it is enough to look at
    # the first factor to find the correct dimension
    for t in product.transforms:
        if t.index0.type == type and t.index0.index == i:
            return (product.factors[0].basisfunction.element.shapedim, True)
        if t.index1.type == type and t.index1.index == i:
            return (product.factors[0].basisfunction.element.shapedim, True)

    # Check Coefficients
    # FIXME: not implemented

    return (0, False)

def dim_factor(factor, i, type):
    """Try to find dimension number i of the given type. If found, the
    dimension is returned as (dim, True). Otherwise (0, False) is
    returned."""

    # Check BasisFunction
    if factor.basisfunction.index.type == type and factor.basisfunction.index.index == i:
        return (factor.basisfunction.element.spacedim, True)

    # Check Derivatives
    for d in factor.derivatives:
        if d.index.type == type and d.index.index == i:
            return (factor.basisfunction.element.shapedim, True)

    # Check that the list of Transforms is empty
    if factor.transforms:
        raise RuntimeError, "Non-empty list of Transforms for Factor."

    return (0, False)

    w1 = u.dx(i)*v.dx(i) + u*v
    w2 = u.dx(i)*v.dx(i) + u*v

def __create_auxiliary(product):
    """Modify secondary indices to auxiliary indices if they don't
    appear inside the integral."""
    p = Product(product)
    for t in p.transforms:
        if not __have_index(product, t.index0):
            t.index0.type = "auxiliary"
        if not __have_index(product, t.index1):
            t.index1.type = "auxiliary"
    # FIXME: Check also coefficients
    return p
    
def __have_index(product, index):
    "Check if the product contains the given Index within a Product."
    for f in product.factors:
        # Check BasisFunction
        i = f.basisfunction.index
        if index.type == i.type and index.index == i.index:
            return True
        # Check Derivatives
        for d in f.derivatives:
            if d.index.type == index.type and d.index.index == index.index:
                return True
        # Check that the list of Transforms is empty
        if f.transforms:
            raise RuntimeError, "Non-empty list of Transforms for Factor."
    return False

if __name__ == "__main__":

    print "Testing reassignment"
    print "--------------------"
    
    element = FiniteElement("Lagrange", 1, "triangle")
    
    u = BasisFunction(element)
    v = BasisFunction(element)
    i = Index()

    w1 = u.dx(i)*v.dx(i) + u*v
    w2 = u.dx(i)*v.dx(i) + u*v

    print w1
    print w2
    
    print reassign_indices(w1)
    print reassign_indices(w2)
