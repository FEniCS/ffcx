__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-10-04 -- 2005-05-05"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# FIAT modules
from FIAT.shapes_new import *
from FIAT.Lagrange import Lagrange, VectorLagrange
from FIAT.DiscontinuousLagrange import DiscontinuousLagrange, DiscontinuousVectorLagrange
from FIAT.CrouzeixRaviart import CrouzeixRaviart
from FIAT.RaviartThomas import RaviartThomas
from FIAT.BDM import BDM
from FIAT.Nedelec import Nedelec

# FFC modules
from dofmap import *
from pointmap import *

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
      shape:  "line", "triangle", "tetrahedron"
      degree: 0, 1, 2, ...

    The shape and degree must match the chosen type of finite element."""

    def __init__(self, name, shape, degree = None, num_components = None):
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
            self.element = Lagrange(fiat_shape, degree)
        elif name == "Vector Lagrange":
            self.element = VectorLagrange(fiat_shape, degree, num_components)
        elif name == "Discontinuous Lagrange":
            print "Warning: element untested"
            self.element = DiscontinuousLagrange(fiat_shape, degree)
        elif name == "Discontinuous vector Lagrange":
            print "Warning: element untested"
            self.element = DiscontinuousVectorLagrange(fiat_shape, degree, num_components)
        elif name == "Crouzeix-Raviart":
            print "Warning: element untested"
            self.element = CrouzeixRaviart(fiat_shape)
        elif name == "Raviart-Thomas":
            print "Warning: element untested"
            self.element = RaviartThomas(fiat_shape, degree)
        elif name == "Brezzi-Douglas-Marini":
            print "Warning: element untested"
            self.element = BDM(fiat_shape, degree)
        elif name == "Nedelec":
            print "Warning: element untested"
            self.element = Nedelec(degree)
        else:
            raise RuntimeError, "Unknown finite element: " + str(name)

        # Create dof map
        self.dofmap = DofMap(fiat_shape, self.element.dual_basis())

        # Create point map
        self.pointmap = PointMap(self.element.dual_basis())
        
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
        return dim[self.shape()]

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

    w1 = P1.basis()[0];
    w2 = P2.basis()[0];

    print w1
    print w2
