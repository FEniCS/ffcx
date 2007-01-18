__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-10-04 -- 2007-01-18"
__copyright__ = "Copyright (C) 2004-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# Modified by Garth N. Wells 2006
# Modified by Marie E. Rognes 2006
# Modified by Andy R Terrel 2007

# Python modules
import sys

# FIAT modules
from FIAT.shapes import *
from FIAT.Lagrange import Lagrange, VectorLagrange
from FIAT.DiscontinuousLagrange import DiscontinuousLagrange, DiscontinuousVectorLagrange
from FIAT.CrouzeixRaviart import CrouzeixRaviart, VectorCrouzeixRaviart
from FIAT.RaviartThomas import RaviartThomas
from FIAT.BDM import BDM
from FIAT.Nedelec import Nedelec

# FFC common modules
sys.path.append("../../")
from ffc.common.debug import *

# FFC compiler modules
from ffc.compiler.nodemap import *
from ffc.compiler.pointmap import *
from ffc.compiler.vertexeval import *

# FFC fem modules
import mixedelement

shape_to_string = {LINE: "line", TRIANGLE: "triangle", TETRAHEDRON: "tetrahedron"}
string_to_shape = {"line": LINE, "triangle": TRIANGLE, "tetrahedron": TETRAHEDRON}

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


    The shape and degree must match the chosen type of finite element.

    The attribute mapping specifies the type of mapping of the
    reference basis to the global basis functions:

      mapping: "Standard", "Piola"

    The transform defaults to "Standard" for H1 and L2 elements,
    "Piola" for H(div) elements.
    """

    # FIXME: Make functions attributes: space_dimension() --> space_dimension etc

    def __init__(self, type, shape, degree = None, num_components = None, element = None):
        "Create FiniteElement."

        # Initialize data
        self.type_str = type
        self.shape_str = shape
        self.element = None
        self.mapping = "Standard" # Default

        if not element is None:
            self.element = element
            self.fiat_shape = element.domain_shape()
        else:
            # Choose shape
            self.fiat_shape = string_to_shape[shape]

            # Choose function space
            if type == "Lagrange":
                self.element = Lagrange(self.fiat_shape, degree)
            elif type == "Vector Lagrange":
                self.element = VectorLagrange(self.fiat_shape, degree, num_components)
            elif type == "Discontinuous Lagrange":
                self.element = DiscontinuousLagrange(self.fiat_shape, degree)
            elif type == "Discontinuous vector Lagrange":
                self.element = DiscontinuousVectorLagrange(self.fiat_shape, degree, num_components)
            elif type == "Crouzeix-Raviart":
                self.element = CrouzeixRaviart(self.fiat_shape)
            elif type == "Vector Crouzeix-Raviart":
                self.element = VectorCrouzeixRaviart(self.fiat_shape)
            elif type == "Raviart-Thomas" or type == "RT":
                print "Warning: element untested"
                self.element = RaviartThomas(self.fiat_shape, degree)
                self.mapping = "Piola"
            elif type == "Brezzi-Douglas-Marini" or type == "BDM":
                print "Warning: element untested"
                self.element = BDM(self.fiat_shape, degree)
                self.mapping = "Piola"
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

        # FIXME: Attributes in new format
        self.fiat_element = self.element
        self.dof_entities = self.fiat_element.dual_basis().entity_ids

        return

    def basis(self):
        "Return basis of finite element space."
        return self.element.function_space()

    def degree(self):
        "Return degree of polynomial basis."
        return self.basis().degree()

    def cell_shape(self):
        "Return the cell shape"
        return self.element.domain_shape()

    def facet_shape(self):
        "Return shape of facet."
        return self.cell_shape() - 1

    def space_dimension(self):
        "Return the dimension of the finite element function space"
        return len(self.basis())

    def shapedim(self):
        "Return dimension of of shape."
        return dim[self.cell_shape()]

    def value_rank(self):
        "Return the rank of the value space"
        return self.basis().rank()
    
    def value_dimension(self, i):
        "Return the dimension of the value space for axis i"
        if self.value_rank() == 0:
            return 1
        else:
            return self.basis().tensor_dim()[i]

    def num_sub_elements(self):
        "Return the number of sub elements"
        return 1

    def sub_element(self, i):
        "Return sub element i"
        return self

    def vectordim(self):
        "Return vector dimension (number of components)"
        if self.value_rank() == 0:
            return 1
        elif self.value_rank() == 1:
            return self.value_dimension(0)
        else:
            raise RuntimeError, "Can only handle scalar or vector-valued elements."

    def num_facets(self):
        "Return number of facets for shape of element."
        if self.element.domain_shape() == TRIANGLE:
            return 3
        elif self.element.domain_shape() == TETRAHEDRON:
            return 4
        else:
            raise RuntimeError, "Unknown shape."

    def num_alignments(self):
        "Return number of possible alignments of two cells at a common facet."
        if self.element.domain_shape() == TRIANGLE:
            return 2
        elif self.element.domain_shape() == TETRAHEDRON:
            return 6
        else:
            raise RuntimeError, "Unknown shape."

    def tabulate(self, order, points, facet = None):
        """Return tabulated values of derivatives up to given order of
        basis functions at given points. If facet is not None, then the
        values are tabulated on the given facet, with the points given
        on the corresponding reference facet."""
        if facet == None:
            return self.element.function_space().tabulate_jet(order, points)
        else:
            facet_shape = self.facet_shape()
            return self.element.function_space().trace_tabulate_jet(facet_shape, facet, order, points)

    def signature(self):
        "Return a string identifying the finite element"
        if self.vectordim() > 1:
            return "%s finite element of degree %d on a %s with %d components" % \
                   (self.type_str, self.degree(), self.shape_str, self.vectordim())
        else:
            return "%s finite element of degree %d on a %s" % \
                   (self.type_str, self.degree(), self.shape_str)

    def __add__(self, other):
        "Create mixed element."
        if isinstance(other, FiniteElement):
            return mixedelement.MixedElement([self, other])
        elif isinstance(other, mixedelement.MixedElement):
            return mixedelement.MixedElement([self] + other.elements)
        else:
            raise RuntimeError, "Unable to create mixed element from given object: " + str(other)

    def __str__(self):
        "Pretty print"
        return self.signature()

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
