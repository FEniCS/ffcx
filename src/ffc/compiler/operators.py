"""This module extends the form algebra with a collection of operators
based on the basic form algebra operations."""

__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2005-09-07 -- 2005-09-08"
__copyright__ = "Copyright (c) 2005 Anders Logg"
__license__  = "GNU GPL Version 2"

# Python modules
import Numeric

# FFC compiler modules
from index import *
from algebra import *

def vec(v):
    "Create vector of scalar functions from given vector-valued function."
    # Check if we already have a vector
    if isinstance(v, list):
        return v
    # Check that function is vector-valued
    if not __rank(v) == 1:
        raise RuntimeError, "Object is not vector-valued."
    # Get vector dimension
    n = __tensordim(v, 0)
    # Create list of scalar components
    return [v[i] for i in range(n)]

def dot(v, w):
    "Return scalar product of v and w."
    # Treat differently, depending on the type of arguments
    if isinstance(v, list) and isinstance(w, list):
        return Numeric.vdot(v, w)
    elif isinstance(v, list):
        return Numeric.vdot(v, vec(w))
    elif isinstance(w, list):
        return Numeric.vdot(vec(v), w)
    # Otherwise, use index-notation
    i = Index()
    return v[i]*w[i]

def grad(v):
    "Return gradient of given function."
    # Get shape dimension
    d = __shapedim(v)
    # Compute gradient
    return [v.dx(i) for i in range(d)]

def div(v):
    "Return divergence of given function."
    # Treat differently, depending on the type of arguments
    if isinstance(v, list):
        # Get shape dimension
        d = __shapedim(v[0])
        # Get vector dimension
        n = len(v)
        # Divergence only defined for n = d
        if not n == d:
            raise RuntimeError, "Divergence only defined for v : R^d --> R^d"
        # Compute divergence
        return Numeric.sum([v[i].dx(i) for i in range(d)])
    # Otherwise, use index-notation
    i = Index()
    return v[i].dx(i)

def rot(v):
    "Return rotation of given function."
    # Make sure that we get a vector
    if not isinstance(v, list):
        return rot(vec(v))
    # Get shape dimension
    d = __shapedim(v[0])
    # Get vector dimension
    n = len(v)
    # Rotation only defined for n = d = 3
    if not (n == 3 and d == 3):
        raise RuntimeError, "Rotation only defined for v : R^3 --> R^3"
    # Compute rotation
    return [v[2].dx(1) - v[1].dx(2), v[0].dx(2) - v[2].dx(0), v[1].dx(0) - v[0].dx(1)]

def curl(v):
    "Alternative name for rot."
    return rot(v)
     
def __shapedim(v):
    "Return shape dimension for given object."
    if isinstance(v, BasisFunction):
        return v.element.shapedim()
    elif isinstance(v, Product):
        return __shapedim(v.basisfunctions[0])
    elif isinstance(v, Sum):
        return __shapedim(v.products[0])
    elif isinstance(v, Function):
        return __shapedim(Sum(v))
    else:
        raise RuntimeError, "Shape dimension is not defined for object: " + str(v)
    return 0

def __rank(v):
    "Return rank for given object."
    if isinstance(v, BasisFunction):
        return v.element.rank() - len(v.component)
    elif isinstance(v, Product):
        return __rank(v.basisfunctions[0])
    elif isinstance(v, Sum):
        return __rank(v.products[0])
    elif isinstance(v, Function):
        return __rank(Sum(v))
    else:
        raise RuntimeError, "Rank is not defined for object: " + str(v)
    return 0

def __tensordim(v, i):
    "Return size of given dimension for given object."
    if i < 0 or i >= __rank(v):
        raise RuntimeError, "Tensor dimension out of range."
    if isinstance(v, BasisFunction):
        return v.element.tensordim(i + len(v.component))
    elif isinstance(v, Product):
        return __tensordim(v.basisfunctions[0], i)
    elif isinstance(v, Sum):
        return __tensordim(v.products[0], i)
    elif isinstance(v, Function):
        return __tensordim(Sum(v), i)
    else:
        raise RuntimeError, "Tensor dimension is not defined for object: " + str(v)
    return 0

if __name__ == "__main__":

    scalar = FiniteElement("Lagrange", "tetrahedron", 2)
    vector = FiniteElement("Vector Lagrange", "tetrahedron", 2)

    v = BasisFunction(scalar)
    u = BasisFunction(scalar)
    w = Function(scalar)

    V = BasisFunction(vector)
    U = BasisFunction(vector)
    W = Function(vector)
    
    i = Index()
    j = Index()

    dx = Integral()

    print dot(grad(v), grad(u))*dx
    print vec(U)
    print dot(U, V)
    print dot(vec(V), vec(U))
    print dot(U, grad(v))
    print div(U)
    print dot(rot(V), rot(U))
    print div(grad(dot(rot(V), U)))*dx
