__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2007-12-10 -- 2009-08-26"
__copyright__ = "Copyright (C) 2007-2008 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Garth N. Wells 2006-2009

# FFC fem modules
from finiteelement import *
from dofrepresentation import *
from quadrature import *
#from mapping import *
#from finiteelement import AFFINE, CONTRAVARIANT_PIOLA, COVARIANT_PIOLA

# FFC common modules
from ffc.common.log import error

# UFL modules
from ufl.classes import Cell, Measure

class QuadratureElement(FiniteElement):
    """Write description of QuadratureElement"""

    def __init__(self, ufl_element, domain=None):
        "Create QuadratureElement"

        # Save UFL element
        self.__ufl_element = ufl_element

        # Handle restrictions
        if domain and isinstance(domain, Cell):
            error("Restriction of QuadratureElement to Cell is not supported because all dofs are internal to the element.")
        elif domain and isinstance(domain, Measure):
            error("Restriction of QuadratureElement to Measure has not been implemented yet.")

        # Save incoming arguments
        self.__cell_shape = string_to_shape[ufl_element.cell().domain()]
        self.__domain = domain

        # Compute number of points per axis from the degree of the element
        self.__num_axis_points = (ufl_element.degree() + 2) / 2

        # Save element family
        self.__family = "Quadrature"

        # Set rank (is rank = 0 for this element?)
        self._rank = 0

        # Save element degree (constant)
        self.__degree = 0

        # Set mapping to AFFINE (not important, I think, for this element)
        self.__mapping = AFFINE

        # Create quadrature (only interested in points)
        points, weights = make_quadrature(self.__cell_shape, self.__num_axis_points)

        # Save number of quadrature points
        self.__num_quad_points = len(points)

        # Create entity IDs, ripped from FIAT/DiscontinuousLagrange.py
        # Used by formdata.py to create the DofMap
        entity_ids = {}
        for d in range( shape_to_dim[ self.__cell_shape ] ):
            entity_ids[d] = {}
            for e in range(num_entities[self.__cell_shape][d]):
                entity_ids[d][e] = []
        entity_ids[ shape_to_dim[ self.__cell_shape ]] = {}
        entity_ids[ shape_to_dim[ self.__cell_shape ]][ 0 ] = range( self.__num_quad_points )
        self.__entity_dofs = [entity_ids]

        # FIXME: We could change calls to element.dual_basis().pts to
        # element.dof_coordinates() which would only support a limited number of
        # elements.
        # Initialise a dummy dual_basis.
        # Used in dofmap.py", line 158, in __compute_dof_coordinates
        self.__dual_basis = [DofRepresentation("Quadrature", [pt]) for pt in points]

    def __repr__(self):
        "Pretty print"
        return self.signature()

    def cell_shape(self):
        "Return the cell shape"
        return self.__cell_shape

    def degree(self):
        "Return degree of polynomial basis"
        return self.__degree

    def domain(self):
        "Return the domain"
        return self.__domain

    def dual_basis(self):
        "Return dummy dual basis of finite element space"
        return self.__dual_basis

    def entity_dofs(self):
        "Return the mapping from entities to dofs"
        return self.__entity_dofs

    def family(self):
        "Return a string indentifying the finite element family"
        return self.__family

    def mapping(self):
        "Return the type of mapping associated with the element."
        return self.__mapping

    def num_axis_points(self):
        "Return the number of quadrature points per axis as specified by user"
        return self.__num_axis_points

    def signature(self):
        "Return a string identifying the QuadratureElement"
        return repr(self.__ufl_element)

    def space_dimension(self):
        "Return the total number of quadrature points"
        return self.__num_quad_points

    def tabulate(self, order, points):
        """Return the identity matrix of size (num_quad_points, num_quad_points),
        in a format that monomialintegration and monomialtabulation understands."""

        # Derivatives are not defined on a QuadratureElement
        # FIXME: currently this check results in a warning (should be RuntimeError)
        # because otherwise some forms fails if QuadratureElement is used in a
        # mixed element e.g.,
        # element = CG + QuadratureElement
        # (v, w) = BasisFunctions(element)
        # grad(w): this is in error and should result in a runtime error
        # grad(v): this should be OK, but currently will raise a warning because
        # derivatives are tabulated for ALL elements in the mixed element.
        # This issue should be fixed in UFL and then we can switch on the
        # RuntimeError again.
        if order:
#            error("Derivatives are not defined on a QuadratureElement")
            print "\n*** WARNING: Derivatives are not defined on a QuadratureElement,"
            print   "             returning values of basisfunction.\n"

        # Check if (the number of ) incoming points are equal to
        # quadrature points... 
        if not len(points) == self.__num_quad_points:
            error("Points must be equal to coordinates of quadrature points")
            
        # Return the identity matrix of size __num_quad_points in a
        # suitable format for monomialintegration.
        values = numpy.identity(self.__num_quad_points, float)
        table = [{(0,)*self.__cell_shape: values}]
        return table

    def value_rank(self):
        "Return the rank of the value space"
        return self._rank

    def ufl_element(self):
        "Return the UFL element"
        return self.__ufl_element

