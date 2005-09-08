"""This module extends the form algebra with a collection of operators
based on the basic form algebra operations."""

__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2005-09-07"
__copyright__ = "Copyright (c) 2005 Anders Logg"
__license__  = "GNU GPL Version 2"

# Python modules
from Numeric import *

# FFC compiler modules
from index import *
from algebra import *

def grad(u):
    "Return gradient of given function."
    # Get shape dimension
    d = __shapedim(u)
    # Compute gradient
    return [u.dx(i) for i in range(d)]

def div(u):
    "Return divergence of given function."
    # Create a new Index
    i = Index()
    # Compute divergence
    return u[i].dx(i)

def rot(u):
    "Return divergence of given function."
    # Get shape dimension
    d = __shapedim(u)
    # Definition depends on the dimension
    raise RuntimeError, "Operator rot() not yet implemented"

def __shapedim(u):
    "Compute shapedim for given object."
    if isinstance(u, BasisFunction):
        return u.element.shapedim()
    elif isinstance(u, Product):
        return u.basisfunctions[0].element.shapedim()
    elif isinstance(u, Sum):
        return u.products[0].basisfunctions[0].element.shapedim()
    elif isinstance(u, Function):
        return u.element.shapedim()
    else:
        raise RuntimeError, "Shape dimension is not defined for object: " + str(u)
    return 0
