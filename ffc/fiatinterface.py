# Copyright (C) 2009-2010 Kristian B. Oelgaard and Anders Logg
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC.  If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Garth N. Wells, 2009.
# Modified by Marie Rognes, 2009-2010.
#
# First added:  2009-03-06
# Last changed: 2011-01-13

# Python modules
from numpy import array

# UFL and FIAT modules
import ufl
import FIAT

# FFC modules
from ffc.log import debug, error
from ffc.quadratureelement import QuadratureElement as FFCQuadratureElement

from ffc.mixedelement import MixedElement
from ffc.restrictedelement import RestrictedElement
from ffc.enrichedelement import EnrichedElement, SpaceOfReals

# Dictionary mapping from domain (cell) to dimension
from ufl.geometry import domain2dim

# Mapping from dimension to number of mesh sub-entities. (In principle,
# ufl.geometry.domain2num_facets contains the same information, but
# with string keys.)
entities_per_dim = {1: [2, 1], 2: [3, 3, 1], 3: [4, 6, 4, 1]}

# Cache for computed elements
_cache = {}

def reference_cell(dim):
    if isinstance(dim, int):
        return FIAT.ufc_simplex(dim)
    else:
        return FIAT.ufc_simplex(domain2dim[dim])

def create_element(ufl_element):

    # Check cache
    if ufl_element in _cache:
        debug("Reusing element from cache")
        return _cache[ufl_element]

    # Create regular FIAT finite element
    if isinstance(ufl_element, ufl.FiniteElement):
        element = _create_fiat_element(ufl_element)

    # Create mixed element (implemented by FFC)
    elif isinstance(ufl_element, ufl.MixedElement):
        elements = _extract_elements(ufl_element)
        element = MixedElement(elements)

    # Create element union (implemented by FFC)
    elif isinstance(ufl_element, ufl.EnrichedElement):
        elements = [create_element(e) for e in ufl_element._elements]
        element = EnrichedElement(elements)

    # Create restricted element(implemented by FFC)
    elif isinstance(ufl_element, ufl.RestrictedElement):
        element = _create_restricted_element(ufl_element)

    else:
        error("Cannot handle this element type: %s" % str(ufl_element))

    # Store in cache
    _cache[ufl_element] = element

    return element

def _create_fiat_element(ufl_element):
    "Create FIAT element corresponding to given finite element."

    family = ufl_element.family()
    degree = ufl_element.degree()

    # Handle the space of the constant
    if family == "Real":
        constant = _create_fiat_element(ufl.FiniteElement("DG", ufl_element.cell(), 0))
        return SpaceOfReals(constant)

    # FIXME: AL: Should this really be here?
    # Handle QuadratureElement
    if family == "Quadrature":
        return FFCQuadratureElement(ufl_element)

    cell = reference_cell(ufl_element.cell().domain())

    # Handle Bubble element as RestrictedElement of P_{k} to interior
    if family == "Bubble":
        V = FIAT.supported_elements["Lagrange"](cell, degree)
        dim = ufl_element.cell().geometric_dimension()
        return RestrictedElement(V, _indices(V, "interior", dim), None)

    # Check if finite element family is supported by FIAT
    if not family in FIAT.supported_elements:
        error("Sorry, finite element of type \"%s\" are not supported by FIAT.", family)

    # Create FIAT finite element
    ElementClass = FIAT.supported_elements[family]
    if degree is None:
        element = ElementClass(cell)
    else:
        element = ElementClass(cell, degree)

    return element

def create_quadrature(shape, num_points):
    """
    Generate quadrature rule (points, weights) for given shape with
    num_points points in each direction.
    """

    # FIXME: KBO: Can this be handled more elegantly?
    if isinstance(shape, int) and shape == 0 or domain2dim[shape] == 0:
        return ([()], array([1.0,]))

    quad_rule = FIAT.make_quadrature(reference_cell(shape), num_points)
    return quad_rule.get_points(), quad_rule.get_weights()

def map_facet_points(points, facet):
    """
    Map points from the e (UFC) reference simplex of dimension d - 1
    to a given facet on the (UFC) reference simplex of dimension d.
    This may be used to transform points tabulated for example on the
    2D reference triangle to points on a given facet of the reference
    tetrahedron.
    """

    # Special case, don't need to map coordinates on vertices
    dim = len(points[0]) + 1
    if dim == 1:
        return [[(0.0,), (1.0,)][facet]]

    # Vertex coordinates
    vertex_coordinates = \
        {1: ((0.,), (1.,)),
         2: ((0., 0.), (1., 0.), (0., 1.)),
         3: ((0., 0., 0.), (1., 0., 0.),(0., 1., 0.), (0., 0., 1))}

    # Facet vertices
    facet_vertices = \
        {2: ((1, 2), (0, 2), (0, 1)),
         3: ((1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2))}

    # Compute coordinates and map
    coordinates = [vertex_coordinates[dim][v] for v in facet_vertices[dim][facet]]
    new_points = []
    for point in points:
        w = (1.0 - sum(point),) + tuple(point)
        x = tuple(sum([w[i]*array(coordinates[i]) for i in range(len(w))]))
        new_points += [x]

    return new_points

def _extract_elements(ufl_element, domain=None):
    "Recursively extract un-nested list of (component) elements."

    elements = []
    if isinstance(ufl_element, ufl.MixedElement):
        for sub_element in ufl_element.sub_elements():
            elements += _extract_elements(sub_element, domain)
        return elements

    # Handle restricted elements since they might be mixed elements too.
    if isinstance(ufl_element, ufl.RestrictedElement):
        base_element = ufl_element.element()
        restriction = ufl_element.domain_restriction()
        return _extract_elements(base_element, restriction)

    if domain:
        ufl_element = ufl.RestrictedElement(ufl_element, domain)

    elements += [create_element(ufl_element)]
    return elements

def _create_restricted_element(ufl_element):
    "Create an FFC representation for an UFL RestrictedElement."

    if not isinstance(ufl_element, ufl.RestrictedElement):
        error("create_restricted_element expects an ufl.RestrictedElement")

    base_element = ufl_element.element()
    domain = ufl_element.domain_restriction()

    # If simple element -> create RestrictedElement from fiat_element
    if isinstance(base_element, ufl.FiniteElement):
        element = _create_fiat_element(base_element)
        return RestrictedElement(element, _indices(element, domain), domain)

    # If restricted mixed element -> convert to mixed restricted element
    if isinstance(base_element, ufl.MixedElement):
        elements = _extract_elements(base_element, domain)
        return MixedElement(elements)

    error("Cannot create restricted element from %s" % str(ufl_element))

def _indices(element, domain, dim=0):
    "Extract basis functions indices that correspond to domain."

    # FIXME: The domain argument in FFC/UFL needs to be re-thought and
    # cleaned-up.

    # If domain is "interior", pick basis functions associated with
    # cell.
    if domain == "interior" and dim:
        return element.entity_dofs()[dim][0]

    # If domain is a ufl.Cell, pick basis functions associated with
    # the topological degree of the domain and of all lower
    # dimensions.
    if isinstance(domain, ufl.Cell):
        dim = domain.topological_dimension()
        entity_dofs = element.entity_dofs()
        indices = []
        for dim in range(domain.topological_dimension() + 1):
            entities = entity_dofs[dim]
            for (entity, index) in entities.iteritems():
                indices += index
        return indices

    # Just extract all indices to make handling in RestrictedElement
    # uniform.
    elif isinstance(domain, ufl.Measure):
        indices = []
        entity_dofs = element.entity_dofs()
        for dim, entities in entity_dofs.items():
            for entity, index in entities.items():
                indices += index
        return indices

    else:
        error("Restriction to domain: %s, is not supported." % repr(domain))

