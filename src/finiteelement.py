__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-10-04"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# FIAT modules
from FIAT.base import shapes
from FIAT.library.Lagrange import Lagrange
from FIAT.library.Hermite import Hermite

class FiniteElement:

    """A FiniteElement represents a finite element in the classical
    sense as defined by Ciarlet. The actual work is done by FIAT and
    this class serves as the interface to FIAT finite elements from
    within FFC.

    A FiniteElement is specified by giving the name and degree of the
    finite element, together with the name of the reference cell:

      name:   "Lagrange", "Hermite", ...
      degree: 0, 1, 2, ...
      shape:  "line", "triangle", "tetrahedron"

    The degree and shape must match the chosen type of finite
    element."""

    def __init__(self, name, degree, shape):
        "Create FiniteElement."

        # Save data
        self.name = name
        self.degree = degree
        self.shape = shape

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
            self.fiat_space = Lagrange(self.fiat_shape, degree)
        elif name == "Hermite":
            self.fiat_space = Hermite(self.fiat_shape, degree)
        else:
            raise RuntimeError, "Unknown space " + str(name)

        # Get the basis
        self.basis = self.fiat_space.getBasis()

        # Set dimensions
        self.spacedim = len(self.basis)
        self.shapedim = shapes.dims[self.fiat_shape]

        return

    def __repr__(self):
        "Print nicely formatted representation of FiniteElement."
        return str(self.name) + " finite element of degree " + \
               str(self.degree) + " on a " + str(self.shape)

if __name__ == "__main__":

    print "Testing finite element"
    print "----------------------"

    P1 = FiniteElement("Lagrange", 1, "triangle")
    P2 = FiniteElement("Lagrange", 2, "triangle")

    print P1
    print P2
