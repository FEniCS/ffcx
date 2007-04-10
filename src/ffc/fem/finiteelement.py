__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-10-04 -- 2007-04-04"
__copyright__ = "Copyright (C) 2004-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# Modified by Garth N. Wells 2006
# Modified by Marie E. Rognes 2006
# Modified by Andy R Terrel 2007

# Python modules
import sys

# FIAT modules
from FIAT.shapes import *
from FIAT.Lagrange import Lagrange
from FIAT.DiscontinuousLagrange import DiscontinuousLagrange
from FIAT.CrouzeixRaviart import CrouzeixRaviart
from FIAT.RaviartThomas import RaviartThomas
from FIAT.BDM import BDM
from FIAT.Nedelec import Nedelec

# FFC common modules
from ffc.common.debug import *

# FFC fem modules
from mapping import *
import mixedelement

# Dictionaries of basic element data
shape_to_string = {LINE: "line", TRIANGLE: "triangle", TETRAHEDRON: "tetrahedron"}
string_to_shape = {"line": LINE, "triangle": TRIANGLE, "tetrahedron": TETRAHEDRON}
shape_to_dim = {LINE: 1, TRIANGLE: 2, TETRAHEDRON: 3}
shape_to_facet = {LINE: None, TRIANGLE: LINE, TETRAHEDRON: TRIANGLE}
shape_to_num_facets = {LINE: 2, TRIANGLE: 3, TETRAHEDRON: 4}

class FiniteElement:
    """A FiniteElement represents a finite element in the classical
    sense as defined by Ciarlet. The actual work is done by FIAT and
    this class serves as the interface to FIAT finite elements from
    within FFC.

    A FiniteElement is specified by giving the family and degree of the
    finite element, together with the name of the reference cell:

      family: "Lagrange", "Hermite", ...
      shape:  "line", "triangle", "tetrahedron"
      degree: 0, 1, 2, ...

    The shape and degree must match the chosen family of finite element.
    """

    def __init__(self, family, shape, degree = None):
        "Create FiniteElement"

        # Save element family
        self.__family = family

        # Get FIAT element from string
        (self.__fiat_element, self.__mapping) = self.__choose_element(family, shape, degree)

        # Get entity dofs from FIAT element
        self.__entity_dofs = [self.__fiat_element.dual_basis().entity_ids]

    def family(self):
        "Return a string indentifying the finite element family"
        return self.__family

    def signature(self):
        "Return a string identifying the finite element"
        return "%s finite element of degree %d on a %s" % \
               (self.__family, self.degree(), shape_to_string[self.cell_shape()])

    def cell_shape(self):
        "Return the cell shape"
        return self.__fiat_element.domain_shape()

    def space_dimension(self):
        "Return the dimension of the finite element function space"
        return len(self.basis())

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

    def degree(self):
        "Return degree of polynomial basis"
        return self.basis().degree()

    def mapping(self):
        "Return the type of mapping associated with the element"
        return self.__mapping

    def cell_dimension(self):
        "Return dimension of shape"
        return shape_to_dim[self.cell_shape()]

    def facet_shape(self):
        "Return shape of facet"
        return shape_to_facet[self.cell_shape()]

    def num_facets(self):
        "Return number of facets for shape of element"
        return shape_to_num_facets[self.cell_shape()]

    def entity_dofs(self):
        "Return the mapping from entities to dofs"
        return self.__entity_dofs

    def basis(self):
        "Return basis of finite element space"
        return self.__fiat_element.function_space()

    def dual_basis(self):
        "Return dual basis of finite element space"
        return self.__fiat_element.dual_basis()

    def tabulate(self, order, points, facet = None):
        """Return tabulated values of derivatives up to given order of
        basis functions at given points. If facet is not None, then the
        values are tabulated on the given facet, with the points given
        on the corresponding reference facet."""
        if facet == None:
            return self.__fiat_element.function_space().tabulate_jet(order, points)
        else:
            facet_shape = self.facet_shape()
            return self.__fiat_element.function_space().trace_tabulate_jet(facet_shape, facet, order, points)

    def __add__(self, other):
        "Create mixed element"
        return mixedelement.MixedElement([self, other])

    def __choose_element(self, family, shape, degree):
        "Choose FIAT finite element from string"

        # Get FIAT shape from string
        fiat_shape = string_to_shape[shape]
    
        # Choose FIAT function space
        if family == "Lagrange" or family == "CG":
            self.__family = "Lagrange"
            return (Lagrange(fiat_shape, degree),
                    Mapping.AFFINE)

        if family == "Discontinuous Lagrange" or family == "DG":
            self.__family = "Discontinuous Lagrange"
            return (DiscontinuousLagrange(fiat_shape, degree),
                    Mapping.AFFINE)

        if family == "Crouzeix-Raviart" or family == "CR":
            self.__family = "Crouzeix-Raviart"
            return (CrouzeixRaviart(fiat_shape),
                    Mapping.AFFINE)

        if family == "Raviart-Thomas" or family == "RT":
            self.__family = "Raviart-Thomas"
            return (RaviartThomas(fiat_shape, degree),
                    Mapping.PIOLA)

        if family == "Brezzi-Douglas-Marini" or family == "BDM":
            self.__family = "Brezzi-Douglas-Marini"
            return (BDM(fiat_shape, degree),
                    Mapping.PIOLA)

        if family == "Nedelec":
            self.__family = "Nedelec"
            return (Nedelec(degree),
                    Mapping.PIOLA) # ?

        # Unknown element
        raise RuntimeError, "Unknown finite element: " + str(family)

    def __repr__(self):
        "Pretty print"
        return self.signature()
