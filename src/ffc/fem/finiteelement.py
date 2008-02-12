__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-10-04 -- 2007-11-26"
__copyright__ = "Copyright (C) 2004-2007 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Garth N. Wells 2006
# Modified by Marie E. Rognes 2008
# Modified by Andy R Terrel 2007

# Python modules
import sys
import numpy

# FIAT modules
from FIAT.shapes import *
from FIAT.transformedspace import *
from FIAT.Lagrange import Lagrange
from FIAT.DiscontinuousLagrange import DiscontinuousLagrange
from FIAT.CrouzeixRaviart import CrouzeixRaviart
from FIAT.RaviartThomas import RaviartThomas
from FIAT.BDM import BDM
from FIAT.BDFM import BDFM
from FIAT.Nedelec import Nedelec

# FFC common modules
from ffc.common.debug import *

# FFC fem modules
from mapping import *
import mixedelement
import referencecell
from  dofrepresentation import *

# Dictionaries of basic element data
shape_to_string = {LINE: "interval", TRIANGLE: "triangle", TETRAHEDRON: "tetrahedron"}
string_to_shape = {"interval": LINE, "triangle": TRIANGLE, "tetrahedron": TETRAHEDRON}
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

        # Get the transformed (according to mapping) function space:
        self.__transformed_space = self.__transformed_function_space()

        # Get entity dofs from FIAT element
        self.__entity_dofs = [self.__fiat_element.dual_basis().entity_ids]

        # Get the dof identifiers from FIAT element
        self.__dual_basis = self.__create_dof_representation(self.__fiat_element.dual_basis().get_dualbasis_types())
        
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
            # return self.basis().tensor_dim()[i]
            # meg: Need tensor_dim in FIAT.transformedspace
            return self.basis().fspace.tensor_dim()[i]

    def num_sub_elements(self):
        "Return the number of sub elements"
        return 1

    def sub_element(self, i):
        "Return sub element i"
        return self

    def degree(self):
        "Return degree of polynomial basis"
        return self.basis().degree()

    def value_mapping(self, component):
        """Return the type of mapping associated with the i'th
        component of the element"""
        return self.__mapping

    def space_mapping(self, i):
        """Return the type of mapping associated with the i'th basis
        function of the element"""
        return self.__mapping
    
    def value_offset(self, component):
        """Given an absolute component (index), return the associated
        subelement and offset of the component""" 
        return (self, 0)

    def space_offset(self, i):
        """Given a basis function number i, return the associated
        subelement and offset""" 
        return (self, 0)
    
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

    def dual_basis(self):
        "Return the representation dual basis of finite element space"
        return self.__dual_basis

    def basis(self):
        "Return basis of finite element space"
        return self.__transformed_space

    def tabulate(self, order, points):
        """Return tabulated values of derivatives up to given order of
        basis functions at given points."""
        return self.basis().tabulate_jet(order, points)

    def basis_elements(self):
        "Returns a list of all basis elements"
        return [self]

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
                    Mapping.CONTRAVARIANT_PIOLA)

        if family == "Brezzi-Douglas-Marini" or family == "BDM":
            self.__family = "Brezzi-Douglas-Marini"
            return (BDM(fiat_shape, degree),
                    Mapping.CONTRAVARIANT_PIOLA)

        if family == "Brezzi-Douglas-Fortin-Marini" or family == "BDFM":
            self.__family = "Brezzi-Douglas-Fortin-Marini"
            return (BDFM(fiat_shape, degree),
                    Mapping.CONTRAVARIANT_PIOLA)

        if family == "Nedelec":
            self.__family = "Nedelec"
            return (Nedelec(fiat_shape, degree),
                    Mapping.COVARIANT_PIOLA)

        # Unknown element
        raise RuntimeError, "Unknown finite element: " + str(family)

    def __transformed_function_space(self):
        """ Transform the function space onto the chosen reference
        cell according to the given mapping of the finite element."""
        function_space = self.__fiat_element.function_space()
        vertices = referencecell.get_vertex_coordinates(self.cell_dimension())
        if self.__mapping == Mapping.AFFINE:
            return AffineTransformedFunctionSpace(function_space, vertices)
        elif self.__mapping == Mapping.CONTRAVARIANT_PIOLA:
            return PiolaTransformedFunctionSpace(function_space, vertices,
                                                 "div")
        elif self.__mapping == Mapping.COVARIANT_PIOLA:
            return PiolaTransformedFunctionSpace(function_space, vertices,
                                                 "curl")
        else:
            raise FormError, "Unknown transform %s" % str(family)

    def __create_dof_representation(self, list_of_fiat_dofs):
        """ Take the FIAT dof representation and convert it to the ffc
        dof representation including transforming the points from one
        reference element onto the other."""

        dofs = []

        # In order to convert from the FIAT geometry to the UFC
        # geometry we need the pushforward and the jacobian of the
        # geometry mapping. Note that the Jacobian returned by FIAT is
        # the one of the mapping from the physical to the reference
        # element. This is the inverse of the nomenclature usually
        # used in this code.
        pushforward = self.basis().pushforward
        fiat_J = self.basis().get_jacobian()
        if self.__mapping == Mapping.CONTRAVARIANT_PIOLA:
            J = 1.0/numpy.linalg.det(fiat_J)*numpy.transpose(fiat_J)
        elif self.__mapping == Mapping.COVARIANT_PIOLA:
            J = numpy.linalg.inv(fiat_J)

        for fiat_dof_repr in list_of_fiat_dofs:
            (name, fiat_pts, fiat_dirs, fiat_weights) = fiat_dof_repr
            directions = []
            weights = []

            # The points on the FIAT reference element are pushed onto
            # the UFC reference element:
            points = [tuple(pushforward(p)) for p in fiat_pts]
            
            # The direction map according to inverse tranpose of the
            # mapping of the element space
            if not fiat_dirs == None:
                directions = [tuple(numpy.dot(J, d)) for d in fiat_dirs]

            # The integrals are not preserved, so we keep the FIAT
            # weights (for now).
            if not fiat_weights == None:
                weights = [w for w in fiat_weights]

            dofs += [DofRepresentation(name, points, directions, weights)]

        return dofs

    def __repr__(self):
        "Pretty print"
        return self.signature()
