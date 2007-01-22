"""This module extends the form algebra with a collection of operators
based on the basic form algebra operations."""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2005-09-07 -- 2007-01-18"
__copyright__ = "Copyright (C) 2005-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# Modified by Ola Skavhaug, 2005
# Modified by Dag Lindbo, 2006
# Modified by Garth N. Wells 2006
# Modified by Kristian Oelgaard 2006

# Python modules
import sys
import numpy

# FFC common modules
sys.path.append("../../")
from ffc.common.exceptions import *

# FFC fem modules
from ffc.fem.finiteelement import *

# FFC compiler modules
from ffc.compiler.projection import *

# FFC formlanguage modules
from ffc.formlanguage.index import *
from ffc.formlanguage.algebra import *

def Identity(n):
    "Return identity matrix of given size."
    # Let numpy handle the identity
    return numpy.identity(n)

def value_rank(v):
    "Return value rank for given object."
    if isinstance(v, BasisFunction):
        return v.element.value_rank() - len(v.component)
    elif isinstance(v, Product):
        return value_rank(v.basisfunctions[0])
    elif isinstance(v, Sum):
        return value_rank(v.products[0])
    elif isinstance(v, Function):
        return value_rank(Sum(v))
    else:
        return numpy.rank(v)
    return 0

def vec(v):
    "Create vector of scalar functions from given vector-valued function."
    # Check if we already have a vector
    if isinstance(v, list):
        return v
    # Check if we have an element of the algebra
    if isinstance(v, Element):
        # Check that we have a vector
        if not value_rank(v) == 1:
            raise FormError, (v, "Cannot create vector from scalar expression.")
        # Get vector dimension
        n = __value_dimension(v, 0)
        # Create list of scalar components
        return [v[i] for i in range(n)]        
    # Let numpy handle the conversion
    if isinstance(v, numpy.ndarray) and len(v.shape) == 1:
        return v.tolist()
    # Unable to find a proper conversion
    raise FormError, (v, "Unable to convert given expression to a vector,")

def dot(v, w):
    "Return scalar product of given functions."
    # Check ranks
    if value_rank(v) == value_rank(w) == 0:
        # Equivalent to standard inner product
        return v*w
    elif value_rank(v) == value_rank(w) == 1:
        # Check dimensions
        if not len(v) == len(w):
            raise FormError, ((v, w), "Dimensions don't match for scalar product.")
        # Use index notation if possible
        if isinstance(v, Element) and isinstance(w, Element):
            i = Index()
            return v[i]*w[i]
        # Otherwise, use numpy.dot
        return numpy.dot(vec(v), vec(w))
    elif value_rank(v) == value_rank(w) == 2:
        
        # Check dimensions
        if not len(v) == len(w):
            raise FormError, ((v, w), "Dimensions don't match for scalar product.")
        # Compute dot product (:) of matrices
        sum = Sum()
        for i in range(len(v)):
            for j in range(len(v)):
                sum = sum + v[i][j]*w[i][j]
        return sum

def cross(v, w):
    "Return cross product of given functions."
    # Check dimensions
    if not len(v) == len(w):
        raise FormError, ((v, w), "Cross product only defined for vectors in R^3.")
    # Compute cross product
    return [v[1]*w[2] - v[2]*w[1], v[2]*w[0] - v[0]*w[2], v[0]*w[1] - v[1]*w[0]]

def trace(v):
    "Return trace of given matrix"
    # Let numpy handle the trace
    return numpy.trace(v)

def transp(v):
    "Return transpose of given matrix."
    # Let numpy handle the transpose."
    return numpy.transpose(v)

def mult(v, w):
    "Compute matrix-matrix product of given matrices."
    # First, convert to numpy.array (safe for both array and list arguments)
    vv = numpy.array(v)
    ww = numpy.array(w)
    if len(vv.shape) == 0 or len(ww.shape) == 0:
        # One argument is a scalar
        return vv*ww
    if len(vv.shape) == len(ww.shape) == 1:
        # Vector times vector
        return numpy.multiply(vv, ww) 
    elif len(vv.shape) == 2 and (len(ww.shape) == 1 or len(ww.shape) == 2):
        # Matvec or matmat product, use matrixmultiply instead
        return numpy.dot(vv, ww)
    else:
        raise FormError, ((v, w), "Dimensions don't match for multiplication.")

def outer(v,w):
    "Return outer product of vector valued functions, p = v'*w"
    # Check that we got a Function
    if not isinstance(v, Function):
        raise FormError, (v, "Outer products are only defined for Functions.")
    if not isinstance(w, Function):
        raise FormError, (w, "Outer products are only defined for Functions.")
    if not len(v) == len(w):
        raise FormError, ((v, w),"Invalid operand dims in outer product")
    
    vv = vec(v)
    ww = vec(w)
    
    return mult(transp([vv]),[ww])
    
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
    if value_rank(v) == 1:
        return [ [D(v[i], j) for j in range(d)] for i in range(len(v)) ]
    # Otherwise assume we have a scalar
    return [D(v, i) for i in range(d)]

def div(v):
    "Return divergence of given function."
    # Use index notation if possible
    if isinstance(v, Element):
        i = Index()
        return v[i].dx(i)
    # Otherwise, compute the sum explicitly
    sum = Sum()
    for i in range(len(v)):
        sum = sum + D(v[i], i)
    return sum

def rot(v):
    "Return rotation of given function."
    # Check dimensions
    if not len(v) == __shapedim(v) == 3:
        raise FormError, (v, "Rotation only defined for v : R^3 --> R^3")
    # Compute rotation
    return [D(v[2], 1) - D(v[1], 2), D(v[0], 2) - D(v[2], 0), D(v[1], 0) - D(v[0], 1)]

def curl(v):
    "Alternative name for rot."
    return rot(v)

def mean(v):
    "Return mean value of given Function (projection onto piecewise constants)."
    # Check that we got a Function
    if not isinstance(v, Function):
        raise FormError, (v, "Mean values are only supported for Functions.")
    # Different projections needed for scalar and vector-valued elements
    element = v.e0
    if element.value_rank() == 0:
        P0 = FiniteElement("Discontinuous Lagrange", element.shape_str, 0)
        pi = Projection(P0)
        return pi(v)
    else:
        P0 = FiniteElement("Discontinuous vector Lagrange", element.shape_str, 0, element.value_dimension(0))
        pi = Projection(P0)
        return pi(v)

def avg(v):
    "Return the average of v across an interior facet."
    if value_rank(v) == 0:
        return 0.5*(v('+') + v('-'))
    else:
        return [0.5*(v[i]('+') + v[i]('-')) for i in range(len(v))]

def jump(v, n):
    "Return the jump of v with respect to the given normal n across an interior facet."
    if value_rank(v) == 0:
        # v is a scalar and n is a vector
        return [v('+')*n[i]('+') + v('-')*n[i]('-') for i in range(len(n))]
    else:
        # v and n are vectors
        sum = Sum();
        for i in range(len(v)):
            sum = sum + v[i]('+')*n[i]('+') + v[i]('-')*n[i]('-')
        return sum

def __shapedim(v):
    "Return shape dimension for given object."
    if isinstance(v, list):
        # Check that all components have the same shape dimension
        for i in range(len(v) - 1):
            if not __shapedim(v[i]) == __shapedim(v[i + 1]):
                raise FormError, (v, "Components have different shape dimensions.")
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
        raise FormError, (v, "Shape dimension is not defined for given expression.")
    return 0

def __value_dimension(v, i):
    "Return size of given dimension for given object."
    if i < 0 or i >= value_rank(v):
        raise FormError, ((v, i), "Tensor dimension out of range.")
    if isinstance(v, BasisFunction):
        return v.element.value_dimension(i + len(v.component))
    elif isinstance(v, Product):
        return __value_dimension(v.basisfunctions[0], i)
    elif isinstance(v, Sum):
        return __value_dimension(v.products[0], i)
    elif isinstance(v, Function):
        return __value_dimension(Sum(v), i)
    else:
        raise FormError, ((v, i), "Tensor dimension is not defined for given expression.")
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
    print trace(mult(Identity(len(V)), grad(V)))
