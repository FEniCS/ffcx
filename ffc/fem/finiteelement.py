__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-10-04"
__copyright__ = "Copyright (C) 2004-2008 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Garth N. Wells 2006-2009
# Modified by Marie E. Rognes 2008
# Modified by Andy R. Terrel 2007
# Modified by Kristian B. Oelgaard 2009
# Last changed: 2009-12-08

# Python modules
import sys
import numpy

# FIAT modules
#from FIAT.shapes import *

from FIAT.transformedspace import *
from FIAT.Lagrange import Lagrange
from FIAT.DiscontinuousLagrange import DiscontinuousLagrange
from FIAT.CrouzeixRaviart import CrouzeixRaviart
from FIAT.RaviartThomas import RaviartThomas
from FIAT.BDM import BDM
from FIAT.BDFM import BDFM
from FIAT.Nedelec import Nedelec
from FIAT.darcystokes import DarcyStokes

# FFC common modules
from ffc.common.log import error, warning

# FFC fem modules
import mixedelement
import referencecell
from dofrepresentation import DofRepresentation

# UFL modules
from ufl.classes import Cell, Measure
from ufl.objects import dc
from ufl.classes import FiniteElementBase

# Dictionaries of basic element data
ufl_domain2fiat_domain = {"vertex": 0, "interval": LINE, "triangle": TRIANGLE, "tetrahedron": TETRAHEDRON}

# Value mappings
AFFINE = "affine"
CONTRAVARIANT_PIOLA = "contravariant Piola"
COVARIANT_PIOLA = "covariant Piola"
mapping_to_int = {AFFINE:0, CONTRAVARIANT_PIOLA:1, COVARIANT_PIOLA:2}

class FiniteElement(FiniteElementBase):
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

    # TODO: KBO: Look at domain argument in create_element, get from ufl_element and check
    def __init__(self, ufl_element, domain=None):
        "Create FiniteElement"

        # Initialise base class
        FiniteElementBase.__init__(self, ufl_element.family(), ufl_element.cell(),
                                   ufl_element.degree(), ufl_element.value_shape())
        # Save string
        self._repr = repr(ufl_element)

        # Save the domain
        self._domain = domain

        # Get FIAT element from string
        (self._fiat_element, self._mapping) = self.__choose_element(ufl_element.family(),
                                                                      ufl_element.cell().domain(), ufl_element.degree())

        if ufl_element.family() not in ("Quadrature", "QE"):
            # Get the transformed (according to mapping) function space:
            self._transformed_space = self.__transformed_function_space()

            # Get entity dofs from FIAT element
            self._entity_dofs = [self._fiat_element.dual_basis().entity_ids]

            # Get the dof identifiers from FIAT element
            self._dual_basis = self.__create_dof_representation(self._fiat_element.dual_basis().get_dualbasis_types())

            # Dofs that have been restricted (it is a subset of self._entity_dofs)
            self._restricted_dofs = []

            # FIXME: This is just a temporary hack to 'support' tensor elements
            self._rank = self.basis().rank()


        # Handle restrictions
        if domain and isinstance(domain, Cell):
            # Some fair warnings and errors
            # Only restrictions to facets are currently supported
            # TODO: This will not handle quadrilateral, hexahedron and the like.
            if domain.domain() != self.cell().facet_domain():
                error("Restriction of FiniteElement to topological entities other than facets is not supported yet.")
            elif self.family() == "Discontinuous Lagrange":
                error("Restriction of Discontinuous Lagrange elements is not supported because all dofs are internal to the element.")
            elif self.family() != "Lagrange":
                warning("Only restriction of Lagrange elements has been tested.")

            # Remove the dofs that live on the interior of the element.
            # TODO: This must of course change if we will support restrictions
            # to other topological entities than facets.
            # Loop dictionaries (there should be only one here?)
            for entity_dict in self._entity_dofs:
                for key, val in entity_dict.items():
                    if key > domain.topological_dimension():
                        for k,v in val.items():
                            val[k] = []
                        continue
                    # Add dofs to the list of dofs that are restricted (still active)
                    for k,v in val.items():
                        self._restricted_dofs.extend(v)

        elif domain and isinstance(domain, Measure):
            # FIXME: Support for restriction to cracks (dc) is only experimental
            if domain == dc:
                warning("Restriciton to Measure dc is only experimental.")
            else:
                error("Restriction of FiniteElement to Measure has not been implemented yet.")

    def basis(self):
        "Return basis of finite element space"
        # Should be safe w.r.t. restrictions, is only used in this module and
        # evaluate_basis and evaluate_basis_derivatives where it is not abused.
        return self._transformed_space

    def component_element(self, component):
        "Return sub element and offset for given component."
        return (self, 0)

    def dual_basis(self):
        "Return the representation dual basis of finite element space"
        # TODO: The restrictions could be handled by __create_dofs_representations,
        # or is it just as elegant to leave it here?
        # Pick out those dofs that are present after restriction.
        # Assuming same numbering as in entity_dofs (as far as I can tell from
        # FIAT, this should be safe fiat/lagrange.py LagrangeDual())
        # FIXME: Experimental support for dc
        if self.domain_restriction() and self.domain_restriction() != dc:
            new_dofs = []
            for d in self._restricted_dofs:
                new_dofs.append(self._dual_basis[d])
            return new_dofs
        return self._dual_basis

    def entity_dofs(self):
        "Return the mapping from entities to dofs"
        return self._entity_dofs

    def extract_elements(self):
        "Extract list of all recursively nested elements."
        return [self]

    def get_coeffs(self):
        "Return the expansion coefficients from FIAT"
        # TODO: Might be able to propagate this to FIAT?
        # FIXME: Experimental support for dc
        if self.domain_restriction() and self.domain_restriction() != dc:
            # Get coefficients and create new table with only the values
            # associated with the restricted dofs.
            coeffs = self.basis().get_coeffs()
            new_coeffs = []
            for dof in self._restricted_dofs:
                new_coeffs.append(coeffs[dof])
            return numpy.array(new_coeffs)
        return self.basis().get_coeffs()

    def mapping(self):
        "Return the type of mapping associated with the element."
        return self._mapping

    def num_sub_elements(self):
        "Return the number of sub elements"
        return 1

    def space_dimension(self):
        "Return the dimension of the finite element function space"
        # FIXME: Experimental support for dc
        if self.domain_restriction() and self.domain_restriction() != dc:
            return len(self._restricted_dofs)
        return len(self.basis())

    def sub_element(self, i):
        "Return sub element i"
        return self

    def tabulate(self, order, points):
        """Return tabulated values of derivatives up to given order of
        basis functions at given points."""
        # TODO: Might be able to propagate this to FIAT?
        # FIXME: Experimental support for dc
        if self.domain_restriction() and self.domain_restriction() != dc:
            # Get basis values and create new table where only the values
            # associated with the restricted dofs are present
            basis_values = self.basis().tabulate_jet(order, points)
            new_basis = []
            for b in basis_values:
                for k,v in b.items():
                    new_vals = []
                    for dof in self._restricted_dofs:
                        new_vals.append(v[dof])
                    b[k] = numpy.array(new_vals)
                new_basis.append(b)
            return new_basis
        return self.basis().tabulate_jet(order, points)

    def value_dimension(self, i):
        "Return the dimension of the value space for axis i"
        if self.value_rank() == 0:
            return 1
        else:
            # return self.basis().tensor_dim()[i]
            # meg: Need tensor_dim in FIAT.transformedspace
            return self.basis().fspace.tensor_dim()[i]

    def value_rank(self):
        "Return the rank of the value space"
        return self.basis().rank()

    def __choose_element(self, family, shape, degree):
        "Choose FIAT finite element from string"

        # Get FIAT shape from string
        fiat_shape = ufl_domain2fiat_domain[shape]

        # Choose FIAT function space
        if family == "Lagrange" or family == "CG":
            return (Lagrange(fiat_shape, degree), AFFINE)

        if family == "Discontinuous Lagrange" or family == "DG":
            return (DiscontinuousLagrange(fiat_shape, degree), AFFINE)

        if family == "Crouzeix-Raviart" or family == "CR":
            return (CrouzeixRaviart(fiat_shape, degree), AFFINE)

        if family == "Raviart-Thomas" or family == "RT":
            return (RaviartThomas(fiat_shape, degree), CONTRAVARIANT_PIOLA)

        if family == "Brezzi-Douglas-Marini" or family == "BDM":
            return (BDM(fiat_shape, degree), CONTRAVARIANT_PIOLA)

        if family == "Brezzi-Douglas-Fortin-Marini" or family == "BDFM":
            return (BDFM(fiat_shape, degree), CONTRAVARIANT_PIOLA)

        if family == "Nedelec 1st kind H(curl)" or family == "N1curl":
            return (Nedelec(fiat_shape, degree), COVARIANT_PIOLA)

        if family == "Darcy-Stokes" or family == "KLMR":
            if not shape == "triangle":
                error("Sorry, Darcy-Stokes element only available on triangles")
            return (DarcyStokes(degree), CONTRAVARIANT_PIOLA)

        if family == "Quadrature" or family == "QE":
            return (None, AFFINE)

        # Unknown element
        error("Unknown finite element: " + str(family))

    def __create_dof_representation(self, list_of_fiat_dofs):
        """ Take the FIAT dof representation and convert it to the ffc
        dof representation including transforming the points from one
        reference element onto the other."""

        # TODO: Modify this for restrictions? See dual_basis(), might be able
        # to propagate this to FIAT
        dofs = []

        # In order to convert from the FIAT geometry to the UFC
        # geometry we need the pushforward and the jacobian of the
        # geometry mapping. Note that the Jacobian returned by FIAT is
        # the one of the mapping from the physical to the reference
        # element. This is the inverse of the nomenclature usually
        # used in this code.
        pushforward = self.basis().pushforward
        fiat_J = self.basis().get_jacobian()
        if self._mapping == CONTRAVARIANT_PIOLA:
            J = 1.0/numpy.linalg.det(fiat_J)*numpy.transpose(fiat_J)
        elif self._mapping == COVARIANT_PIOLA:
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

    def __transformed_function_space(self):
        """ Transform the function space onto the chosen reference
        cell according to the given mapping of the finite element."""
        function_space = self._fiat_element.function_space()
        vertices = referencecell.get_vertex_coordinates(self.cell().domain())
        if self._mapping == AFFINE:
            return AffineTransformedFunctionSpace(function_space, vertices)
        elif self._mapping == CONTRAVARIANT_PIOLA:
            return PiolaTransformedFunctionSpace(function_space, vertices, "div")
        elif self._mapping == COVARIANT_PIOLA:
            return PiolaTransformedFunctionSpace(function_space, vertices, "curl")
        else:
            error(family, "Unknown transform")

