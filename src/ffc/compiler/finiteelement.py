__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-10-04"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# FIAT modules
from FIAT.base import shapes
from FIAT.library.Lagrange import Lagrange
from FIAT.library.Hermite import Hermite

# FFC modules
from tensorspace import *

class FiniteElement:

    """A FiniteElement represents a finite element in the classical
    sense as defined by Ciarlet. The actual work is done by FIAT and
    this class serves as the interface to FIAT finite elements from
    within FFC. Maybe this class should be moved to FIAT?

    A FiniteElement is specified by giving the name and degree of the
    finite element, together with the name of the reference cell:

      name:   "Lagrange", "Hermite", ...
      degree: 0, 1, 2, ...
      shape:  "line", "triangle", "tetrahedron"
      dims:   [int, int, ...] (optional)

    The degree and shape must match the chosen type of finite element.

    The list dims should contain a list of integers specifying the
    tensor-dimensions of the finite element. The default is an empty
    list, corresponding to a scalar-valued element."""

    def __init__(self, name, shape, degree, dims = []):
        "Create FiniteElement."

        # Allow scalar input of dimension for a vector-valued element
        if isinstance(dims, int): dims = [dims]

        # Initialize data
        self.name = name
        self.shape = shape
        self.degree = degree
        self.fiat_shape = None

        self.space = None
        self.basis = None
        
        self.spacedim = None
        self.shapedim = None
        self.tensordims = None
        self.rank = None
        
        # Choose shape
        if shape == "line":
            self.fiat_shape = shapes.LINE
        elif shape == "triangle":
            self.fiat_shape = shapes.TRIANGLE
        elif shape == "tetrahedron":
            self.fiat_shape = shapes.TETRAHEDRON
        else:
            raise RuntimeError, "Unknown shape " + str(shape)

        # Choose function space
        if name == "Lagrange":
            self.space = Lagrange(self.fiat_shape, degree)
            if dims:
                self.space = TensorSpace(self.space, dims)
        elif name == "Hermite":
            self.space = Hermite(self.fiat_shape, degree)
            if dims:
                self.space = TensorSpace(self.space, dims)
        else:
            raise RuntimeError, "Unknown space " + str(name)

        # Get the basis
        self.basis = self.space.getBasis()

        # Set dimensions
        self.spacedim = len(self.basis)
        self.shapedim = shapes.dims[self.fiat_shape]
        self.tensordims = dims
        self.rank = len(dims)

        return

    def __repr__(self):
        "Print nicely formatted representation of FiniteElement."
        return "%s finite element of degree %d on a %s (rank %d)" % \
               (self.name, self.degree, self.shape, self.rank)

if __name__ == "__main__":

    print "Testing finite element"
    print "----------------------"

    P1 = FiniteElement("Lagrange", "triangle", 1)
    P2 = FiniteElement("Lagrange", "triangle", 2)
    Q5 = FiniteElement("Lagrange", "triangle", 5, 3)

    print P1
    print P2
    print Q5
