"""This module provides a set of FFC built-ins such as special
functions and predefined variables."""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2006-12-01 -- 2007-05-10"
__copyright__ = "Copyright (C) 2006-2007 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# FFC fem modules
from ffc.fem.vectorelement import VectorElement
from ffc.fem.finiteelement import FiniteElement, shape_to_dim, string_to_shape

# FFC language modules
from algebra import Function
from integral import Integral
from index import Index

def Constant(shape):
    """This like a class but is really a function. It returns a DG(0)
    function which may be thought of as a Constant."""

    # Create discontinuous Lagrange element
    element = FiniteElement("Discontinuous Lagrange", shape, 0)
    return Function(element)('+/-')

def VectorConstant(shape, vector_dim=None):
    """This looks like a class but is really a function. It returns a
    vector DG(0) function which may be thought of as a vector
    Constant."""

    # Find out vector dimension
    if vector_dim == None:
        vector_dim = shape_to_dim[string_to_shape[shape]]

    # Create discontinuous vector Lagrange element
    element = VectorElement("Discontinuous Lagrange", shape, 0, vector_dim)
    return Function(element)('+/-')

def MeshSize(shape):
    """This looks like a class but is really a function. It returns a DG(0)
    function which may function as the placeholder for a function defining
    the mesh size."""
    element = FiniteElement("Discontinuous Lagrange", shape, 0)
    h = Function(element)
    h.name = "mesh size"
    return h

def FacetNormal(shape):
    """This looks like a class but is really a function. It returns a
    vectorDG(0) function which may function as the placeholder for a
    function defining the facet normal."""
    if shape == "interval":
        element = FiniteElement("Discontinuous Lagrange", shape, 0)
    else:
        element = VectorElement("Discontinuous Lagrange", shape, 0)
    n = Function(element)
    n.name = "facet normal"
    return n

# Predefined integrals
dx = Integral("cell")
ds = Integral("exterior facet")
dS = Integral("interior facet")

# Predefined indices
i = Index()
j = Index()
k = Index()
l = Index()
m = Index()
n = Index()
