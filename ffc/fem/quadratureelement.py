__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2007-12-10"
__copyright__ = "Copyright (C) 2007-2008 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Garth N. Wells 2006-2009
# Last changed: 2009-12-08

from FIAT.shapes import num_entities

# FFC fem modules
from finiteelement import ufl_domain2fiat_domain
from dofrepresentation import DofRepresentation
from quadrature import make_quadrature
from finiteelement import FiniteElement, AFFINE, CONTRAVARIANT_PIOLA, COVARIANT_PIOLA

# FFC common modules
from ffc.common.log import error

# UFL modules
from ufl.classes import Cell, Measure

import numpy

# Default quadrature element degree
default_quadrature_degree = 1

class QuadratureElement(FiniteElement):
    """Write description of QuadratureElement"""

    def __init__(self, ufl_element, domain=None):
        "Create QuadratureElement"

        # Handle restrictions
        if domain and isinstance(domain, Cell):
            error("Restriction of QuadratureElement to Cell is not supported because all dofs are internal to the element.")
        elif domain and isinstance(domain, Measure):
            error("Restriction of QuadratureElement to Measure has not been implemented yet.")

        FiniteElement.__init__(self, ufl_element, domain)

        # Compute number of points per axis from the degree of the element
        degree = ufl_element.degree()
        if degree is None:
            degree = default_quadrature_degree
        self._num_axis_points = (degree + 2) / 2

        # Set rank (is rank = 0 for this element?)
        self._rank = 0

        # Set element degree to constant
        # FIXME: Is this necessary? At this point we should already have computed
        # the total degree of the form.
        self._degree = 0

        # FIXME: Do we really need to do this here? We just use the number of
        # FIXME: points while the actual points are computed in quadraturerepresentation.py

        # Create quadrature (only interested in points)
        # TODO: KBO: What should we do about quadrature functions that live on ds, dS?
        points, weights = make_quadrature(self.cell().domain(), self._num_axis_points)

        # Save number of quadrature points
        self._num_quad_points = len(points)

        # Create entity IDs, ripped from FIAT/DiscontinuousLagrange.py
        # Used by formdata.py to create the DofMap
        entity_ids = {}
        for d in range( self.cell().topological_dimension() ):
            entity_ids[d] = {}
            for e in range(num_entities[ufl_domain2fiat_domain[self.cell().domain()]][d]):
                entity_ids[d][e] = []
        entity_ids[ self.cell().topological_dimension() ] = {}
        entity_ids[ self.cell().topological_dimension() ][ 0 ] = range( self._num_quad_points )
        self._entity_dofs = [entity_ids]

        # FIXME: We could change calls to element.dual_basis().pts to
        # element.dof_coordinates() which would only support a limited number of
        # elements.
        # Initialise a dummy dual_basis.
        # Used in dofmap.py", line 158, in __compute_dof_coordinates
        self._dual_basis = [DofRepresentation("Quadrature", [pt]) for pt in points]

    def num_axis_points(self):
        "Return the number of quadrature points per axis as specified by user"
        return self._num_axis_points

    def space_dimension(self):
        "Return the total number of quadrature points"
        return self._num_quad_points

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

        # Check that incoming points are as many as the quadrature points
        if not len(points) == self._num_quad_points:
            error("Points must be equal to coordinates of quadrature points")

        # Return the identity matrix of size __num_quad_points in a
        # suitable format for monomialintegration.
        values = numpy.identity(self._num_quad_points, float)
        table = [{(0,)*self.cell().topological_dimension(): values}]
        return table

    def value_rank(self):
        "Return the rank of the value space"
        return self._rank

