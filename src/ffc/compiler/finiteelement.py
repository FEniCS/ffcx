__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-10-04"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# FIAT modules
from FIAT import shapes
from FIAT import quadrature
from FIAT.Lagrange import Lagrange
#from FIAT import Hermite

# FFC compiler  modules
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

        self.basis = None
        self.table = None
        
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
            self.basis = Lagrange(self.fiat_shape, degree)
            if dims:
                tensorspace = TensorSpace(self.basis, dims)
                self.basis = tensorspace.basis
        elif name == "Hermite":
            raise RuntimeError, "Hermite elements not implemented by FIAT."
            self.basis = Hermite(self.fiat_shape, degree)
            if dims:
                tensorspace = TensorSpace(self.basis, dims)
                self.basis = tensorspace.basis
        else:
            raise RuntimeError, "Unknown space " + str(name)

        # Set dimensions
        self.spacedim = len(self.basis)
        self.shapedim = shapes.dims[self.fiat_shape]
        self.tensordims = dims
        self.rank = len(dims)

        return

    # FIXME: Not used, should be moved to FIAT
    def tabulate(self, points, dmax):
        """Tabulate the values of all basis functions and their
        derivatives up to order dmax at all quadrature points."""

        if self.table:
            print "Already tabulated, skipping"

        if dmax > 1:
            debug("Warning: only tabulating first derivatives so this might be slow.")

        # Tabulate basis functions
        self.table = [self.basis.tabulate(points)]

        # Tabulate first derivatives of basis functions
        if self.fiat_shape == shapes.TRIANGLE:
            self.table += [self.basis.deriv_all(0).tabulate(points), \
                           self.basis.deriv_all(1).tabulate(points)]
        elif self.fiat_shape == shapes.TETRAHEDRON:
            self.table += [self.basis.deriv_all(0).tabulate(points), \
                           self.basis.deriv_all(1).tabulate(points), \
                           self.basis.deriv_all(2).tabulate(points)]
        else:
            raise RuntimeError, "Unsupported shape."

    def __repr__(self):
        "Print nicely formatted representation of FiniteElement."
        return "%s of degree %d on a %s (rank %d)" % \
               (self.name, self.degree, self.shape, self.rank)

if __name__ == "__main__":

    print "Testing finite element"
    print "----------------------"

    P1 = FiniteElement("Lagrange", "triangle", 1)
    Q1 = FiniteElement("Lagrange", "triangle", 1, 3)
    
    P2 = FiniteElement("Lagrange", "triangle", 2)

    print P1
    print Q1
    print P2

    quadrature = quadrature.make_quadrature(P1.fiat_shape, 5)

    w1 = P1.basis[0];
    w2 = Q1.basis[0][0];

    I1 = quadrature(w1.deriv(0))
    I2 = quadrature(w2.deriv(0))

    print "I1 = " + str(I1)
    print "I2 = " + str(I2)
