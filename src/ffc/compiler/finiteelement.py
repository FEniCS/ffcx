__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-10-04 -- 2006-02-20"
__copyright__ = "Copyright (C) 2004-2006 Anders Logg"
__license__  = "GNU GPL Version 2"

# Python modules
import sys

# FIAT modules
from FIAT.shapes import *
from FIAT.Lagrange import Lagrange, VectorLagrange
from FIAT.DiscontinuousLagrange import DiscontinuousLagrange, DiscontinuousVectorLagrange
from FIAT.CrouzeixRaviart import CrouzeixRaviart
from FIAT.RaviartThomas import RaviartThomas
from FIAT.BDM import BDM
from FIAT.Nedelec import Nedelec

# FFC common modules
sys.path.append("../../")
from ffc.common.debug import *

# FFC compiler modules
from nodemap import *
from pointmap import *
from vertexeval import *
import mixedelement

shape_to_string = { LINE: "Line", TRIANGLE: "Triangle", TETRAHEDRON: "Tetrahedron" }
string_to_shape = { "Line": LINE, "Triangle": TRIANGLE, "Tetrahedron": TETRAHEDRON }

class FiniteElement:

    """A FiniteElement represents a finite element in the classical
    sense as defined by Ciarlet. The actual work is done by FIAT and
    this class serves as the interface to FIAT finite elements from
    within FFC.

    A FiniteElement is specified by giving the type and degree of the
    finite element, together with the name of the reference cell:

      type:   "Lagrange", "Hermite", ...
      shape:  "line", "triangle", "tetrahedron"
      degree: 0, 1, 2, ...

    The shape and degree must match the chosen type of finite element."""

    def __init__(self, type, shape, degree = None, num_components = None, element = None):
        "Create FiniteElement."

        # Initialize data
        self.type_str = type
        self.shape_str = shape
        self.element = None

        if not element is None:
            self.element = element
            self.fiat_shape = element.domain_shape()
        else:
            # Choose shape
            if shape == "line":
                self.fiat_shape = LINE
            elif shape == "triangle":
                self.fiat_shape = TRIANGLE
            elif shape == "tetrahedron":
                self.fiat_shape = TETRAHEDRON
            else:
                raise RuntimeError, "Unknown shape " + str(shape)

            # Choose function space
            if type == "Lagrange":
                self.element = Lagrange(self.fiat_shape, degree)
            elif type == "Vector Lagrange":
                self.element = VectorLagrange(self.fiat_shape, degree, num_components)
            elif type == "Discontinuous Lagrange":
                self.element = DiscontinuousLagrange(self.fiat_shape, degree)
            elif type == "Discontinuous vector Lagrange":
                self.element = DiscontinuousVectorLagrange(self.fiat_shape, degree, num_components)
                # Check for known FIAT bug
                if (not num_components == None) and (not num_components == self.element.function_space().tensor_dim()[0]):
                    raise RuntimeError, \
"""Discontinous vector Lagrange element has wrong number of components.
You need to patch your installation of FIAT. For more information, see
http://www.fenics.org/pipermail/fiat-dev/2005-August/000060.html"""
            elif type == "Crouzeix-Raviart":
                print "Warning: element untested"
                self.element = CrouzeixRaviart(self.fiat_shape)
            elif type == "Raviart-Thomas":
                print "Warning: element untested"
                self.element = RaviartThomas(self.fiat_shape, degree)
            elif type == "Brezzi-Douglas-Marini":
                print "Warning: element untested"
                self.element = BDM(self.fiat_shape, degree)
            elif type == "Nedelec":
                print "Warning: element untested"
                self.element = Nedelec(degree)
            else:
                raise RuntimeError, "Unknown finite element: " + str(type)

        # Save dual basis
        self.fiat_dual = self.element.dual_basis()

        # Create node map
        self.nodemap = NodeMap(self)

        # Create point map
        self.pointmap = PointMap(self)

        # Create vertex evaluation
        self.vertexeval = VertexEval(self)

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

    def vectordim(self):
        "Return vector dimension (number of components)"
        if self.rank() == 0:
            return 1
        elif self.rank() == 1:
            return self.tensordim(0)
        else:
            raise RuntimeError, "Can only handle scalar or vector-valued elements."

    def tabulate(self, order, points):
        """Return tabulated values of derivatives up to given order of
        basis functions at given points."""
        return self.element.function_space().tabulate_jet(order, points)

    def __add__(self, other):
        "Create mixed element."
        if isinstance(other, FiniteElement):
            return mixedelement.MixedElement([self, other])
        elif isinstance(other, mixedelement.MixedElement):
            return mixedelement.MixedElement([self] + other.elements)
        else:
            raise RuntimeError, "Unable to create mixed element from given object: " + str(other)

    def __repr__(self):
        "Print nicely formatted representation of FiniteElement."
        if self.vectordim() > 1:
            return "%s finite element of degree %d on a %s with %d components" % \
                   (self.type_str, self.degree(), self.shape_str, self.vectordim())
        else:
            return "%s finite element of degree %d on a %s" % \
                   (self.type_str, self.degree(), self.shape_str)

if __name__ == "__main__":

    print "Testing finite element"
    print "----------------------"

    P1 = FiniteElement("Lagrange", "triangle", 1)
    P2 = FiniteElement("Lagrange", "triangle", 2)

    w1 = P1.basis()[0];
    w2 = P2.basis()[0];

    x0 = (-1.0, -1.0)
    x1 = (1.0, -1.0)
    x2 = (-1.0, 1.0)

    x3 = (-0.79456469038074751, -0.82282408097459203)

    print w1(x0), w1(x1), w1(x2)
    print w2(x0), w2(x1), w2(x2)

    print w1(x3)
    print w2(x3)
