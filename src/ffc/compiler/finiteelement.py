__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-10-04 -- 2005-05-02"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# FIAT modules
from FIAT import quadrature
from FIAT.shapes import *
from FIAT.Lagrange import Lagrange, VectorLagrange, NewVectorLagrange

# FFC modules
from dofmap import *

shape_to_string = { LINE: "Line", TRIANGLE: "Triangle", TETRAHEDRON: "Tetrahedron" }
string_to_shape = { "Line": LINE, "Triangle": TRIANGLE, "Tetrahedron": TETRAHEDRON }

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

    The degree and shape must match the chosen type of finite element."""

    def __init__(self, name, shape, degree, num_components = None):
        "Create FiniteElement."

        # Initialize data
        self.name = name
        self.element = None

        # Choose shape
        if shape == "line":
            fiat_shape = LINE
        elif shape == "triangle":
            fiat_shape = TRIANGLE
        elif shape == "tetrahedron":
            fiat_shape = TETRAHEDRON
        else:
            raise RuntimeError, "Unknown shape " + str(shape)

        # Choose function space
        if name == "Lagrange":

            # Sanity check
            if num_components:
                raise RuntimeError("Number of components cannot be specified for scalar element (use \"Vector Lagrange\")")

        # Choose function space

            self.element = Lagrange(fiat_shape, degree)

        elif name == "Vector Lagrange":

            self.element = NewVectorLagrange(fiat_shape, degree, num_components)

        else:
            raise RuntimeError, "Unknown finite element: " + str(name)

        # Create dof map
        self.dofmap = DofMap(fiat_shape, self.element.dual_basis())

        return

    def basis(self):
        "Return basis of finite element space."
        return self.element.function_space()

    def degree(self):
        "Return degree of polynomial basis."
        return self.basis().degree()

    def shape(self):
        "Return shape used for element."
        return self.element.domain_shape()

    def spacedim(self):
        "Return dimension of finite element space."
        return len(self.basis())

    def shapedim(self):
        "Return dimension of of shape."
        return dims[self.shape()]

    def rank(self):
        "Return rank of basis functions."
        return self.basis().rank()
    
    def tensordim(self, i):
        "Return size of given dimension."
        tensordims = self.basis().tensor_dim()
        return tensordims[i]

    def __repr__(self):
        "Print nicely formatted representation of FiniteElement."
        return "%s finite element of degree %d on a %s" % \
               (self.name, self.degree(), shape_to_string[self.shape()])

if __name__ == "__main__":

    print "Testing finite element"
    print "----------------------"

    P1 = FiniteElement("Lagrange", "triangle", 1)
    P2 = FiniteElement("Lagrange", "triangle", 2)

    quadrature = quadrature.make_quadrature(P1.shape(), 5)

    w1 = P1.basis()[0];
    w2 = P2.basis()[0];

    I1 = quadrature(w1.deriv(0))
    I2 = quadrature(w2.deriv(0))

    print "I1 = " + str(I1)
    print "I2 = " + str(I2)
