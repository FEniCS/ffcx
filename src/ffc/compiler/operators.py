"""This module extends the form algebra with a collection of operators
based on the basic form algebra operations."""

__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2005-09-07 -- 2005-09-09"
__copyright__ = "Copyright (c) 2005 Anders Logg"
__license__  = "GNU GPL Version 2"

# Python modules
import sys
import Numeric

# FFC compiler modules
sys.path.append("../../")
from index import *
from algebra import *

def Identity(n):
    "Return identity matrix of given size."
    # Let Numeric handle the identity
    return Numeric.identity(n)

def rank(v):
    "Return rank for given object."
    if isinstance(v, BasisFunction):
        return v.element.rank() - len(v.component)
    elif isinstance(v, Product):
        return rank(v.basisfunctions[0])
    elif isinstance(v, Sum):
        return rank(v.products[0])
    elif isinstance(v, Function):
        return rank(Sum(v))
    else:
        return Numeric.rank(v)
    return 0

def vec(v):
    "Create vector of scalar functions from given vector-valued function."
    # Check if we already have a vector
    if isinstance(v, list):
        return v
    # Check that function is vector-valued
    if not rank(v) == 1:
        raise RuntimeError, "Object is not vector-valued."
    # Get vector dimension
    n = __tensordim(v, 0)
    # Create list of scalar components
    return [v[i] for i in range(n)]

def dot(v, w):
    "Return scalar product of given functions."
    # Check ranks
    if rank(v) == rank(w) == 1:
        # Check dimensions
        if not len(v) == len(w):
            raise RuntimeError, "Dimensions don't match for scalar product."
        # Use index notation if possible
        if isinstance(v, Element) and isinstance(w, Element):
            i = Index()
            return v[i]*w[i]
        # Otherwise, use Numeric.vdot
        return Numeric.vdot(vec(v), vec(w))
    elif rank(v) == rank(w) == 2:
        # Check dimensions
        if not len(v) == len(w):
            raise RuntimeError, "Dimensions don't match for scalar product."
        # Compute dot product (:) of matrices
        return Numeric.sum([v[i][j]*w[i][j] for i in range(len(v)) for j in range(len(v[i]))])

def cross(v, w):
    "Return cross product of given functions."
    # Check dimensions
    if not len(v) == len(w):
        raise RuntimeError, "Cross product only defined for vectors in R^3."
    # Compute cross product
    return [v[1]*w[2] - v[2]*w[1], v[2]*w[0] - v[0]*w[2], v[0]*w[1] - v[1]*w[0]]

def trace(v):
    "Return trace of given matrix"
    # Let Numeric handle the trace
    return Numeric.trace(v)

def transp(v):
    "Return transpose of given matrix."
    # Let Numeric handle the transpose."
    return Numeric.transpose(v)

def mult(v, w):
    "Compute matrix-matrix product of given matrices."
    # Let Numeric handle the product."
    return Numeric.multiply(v, w)

def D(v, i):
    "Return derivative of v in given coordinate direction."
    # Use member function dx() if possible
    if isinstance(v, Element):
        return v.dx(i)
    # Otherwise, apply to each component
    return [D(v[j], i) for j in range(len(v))]
    
def grad(v):
    "Return gradient of given function."
    # Get shape dimension
    d = __shapedim(v)
    # Check if we have a vector
    if rank(v) == 1:
        return [ [D(v[i], j) for j in range(d)] for i in range(len(v)) ]
    # Otherwise assume we have a scalar
    return [D(v, i) for i in range(d)]

def div(v):
    "Return divergence of given function."
    # Check dimensions
    if not len(v) == __shapedim(v):
        raise RuntimeError, "Dimensions don't match for divergence."
    # Use index notation if possible
    if isinstance(v, Element):
        i = Index()
        return v[i].dx(i)
    # Otherwise, use Numeric.sum
    return Numeric.sum([D(v[i], i) for i in range(len(v))])

def rot(v):
    "Return rotation of given function."
    # Check dimensions
    if not len(v) == __shapedim(v) == 3:
        raise RuntimeError, "Rotation only defined for v : R^3 --> R^3"
    # Compute rotation
    return [D(v[2], 1) - D(v[1], 2), D(v[0], 2) - D(v[2], 0), D(v[1], 0) - D(v[0], 1)]

def curl(v):
    "Alternative name for rot."
    return rot(v)
     
def __shapedim(v):
    "Return shape dimension for given object."
    if isinstance(v, list):
        # Check that all components have the same shape dimension
        for i in range(len(v) - 1):
            if not __shapedim(v[i]) == __shapedim(v[i + 1]):
                raise RuntimeError, "Components have different shape dimensions."
        # Return length of first term
        return __shapedim(v[0])
    elif isinstance(v, BasisFunction):
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

def __tensordim(v, i):
    "Return size of given dimension for given object."
    if i < 0 or i >= rank(v):
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

    i = Index()

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
    print cross(V, U)
    print trace(mult(I(len(V)), grad(V)))
