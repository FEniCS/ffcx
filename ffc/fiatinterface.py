# Copyright (C) 2009-2013 Kristian B. Oelgaard and Anders Logg
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Garth N. Wells, 2009.
# Modified by Marie Rognes, 2009-2013.
# Modified by Martin Alnaes, 2013
#
# First added:  2009-03-06
# Last changed: 2013-01-25

# Python modules
from numpy import array

# UFL and FIAT modules
import ufl
import FIAT

# FFC modules
from ffc.log import debug, error, ffc_assert
from ffc.quadratureelement import QuadratureElement as FFCQuadratureElement
from ffc.timeelements import LobattoElement as FFCLobattoElement
from ffc.timeelements import RadauElement as FFCRadauElement

from ffc.mixedelement import MixedElement
from ffc.restrictedelement import RestrictedElement
from ffc.enrichedelement import EnrichedElement, SpaceOfReals

# Dictionary mapping from cellname to dimension
from ufl.geometry import cellname2dim

# Number of entities associated with each cell name
cellname_to_num_entities = {
    "cell1D": None,
    "cell2D": None,
    "cell3D": None,
    "interval": (2, 1),
    "triangle": (3, 3, 1),
    "tetrahedron": (4, 6, 4, 1),
    "quadrilateral": (4, 4, 1),
    "hexahedron": (8, 12, 6, 1),
    }

# Element families supported by FFC
supported_families = ("Brezzi-Douglas-Marini",
                      "Brezzi-Douglas-Fortin-Marini",
                      "Crouzeix-Raviart",
                      "Discontinuous Lagrange",
                      "Lagrange",
                      "Lobatto",
                      "Nedelec 1st kind H(curl)",
                      "Nedelec 2nd kind H(curl)",
                      "Radau",
                      "Raviart-Thomas",
                      "Real",
                      "Bubble",
                      "Quadrature")

# Mapping from dimension to number of mesh sub-entities. (In principle,
# cellname_to_num_entities contains the same information, but with string keys.)
# DISABLED ON PURPOSE: It's better to use cell name instead of dimension
#                      to stay generic w.r.t. future box elements.
#entities_per_dim = {1: [2, 1], 2: [3, 3, 1], 3: [4, 6, 4, 1]}

# Cache for computed elements
_cache = {}

def reference_cell(dim):
    if isinstance(dim, int):
        return FIAT.ufc_simplex(dim)
    else:
        return FIAT.ufc_simplex(cellname2dim[dim])

def reference_cell_vertices(dim):
    "Return dict of coordinates of reference cell vertices for this 'dim'."
    cell = reference_cell(dim)
    return cell.get_vertices()

def create_element(ufl_element):

    # Create element signature for caching (just use UFL element)
    element_signature = ufl_element

    # Check cache
    if element_signature in _cache:
        debug("Reusing element from cache")
        return _cache[element_signature]

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
    _cache[element_signature] = element

    return element

def _create_fiat_element(ufl_element):
    "Create FIAT element corresponding to given finite element."

    # Get element data
    family = ufl_element.family()
    cell = ufl_element.cell()
    degree = ufl_element.degree()

    # Check that FFC supports this element
    ffc_assert(family in supported_families,
               "This element family (%s) is not supported by FFC." % family)

    # Handle the space of the constant
    if family == "Real":
        dg0_element = ufl.FiniteElement("DG", cell, 0)
        constant = _create_fiat_element(dg0_element)
        element = SpaceOfReals(constant)

    # Handle the specialized time elements
    elif family == "Lobatto" :
        element = FFCLobattoElement(ufl_element.degree())

    elif family == "Radau" :
        element = FFCRadauElement(ufl_element.degree())

    # FIXME: AL: Should this really be here?
    # Handle QuadratureElement
    elif family == "Quadrature":
        element = FFCQuadratureElement(ufl_element)

    else:
        # Create FIAT cell
        fiat_cell = reference_cell(cell.cellname())

        # Handle Bubble element as RestrictedElement of P_{k} to interior
        if family == "Bubble":
            V = FIAT.supported_elements["Lagrange"](fiat_cell, degree)
            dim = cell.topological_dimension()
            return RestrictedElement(V, _indices(V, "interior", dim), None)

        # Check if finite element family is supported by FIAT
        if not family in FIAT.supported_elements:
            error("Sorry, finite element of type \"%s\" are not supported by FIAT.", family)

        # Create FIAT finite element
        ElementClass = FIAT.supported_elements[family]
        if degree is None:
            element = ElementClass(fiat_cell)
        else:
            element = ElementClass(fiat_cell, degree)

    # Consistency check between UFL and FIAT elements. This will not hold for elements
    # where the reference value shape is different from the global value shape, i.e.
    # RT elements on a triangle in 3D.
    #ffc_assert(element.value_shape() == ufl_element.value_shape(),
    #           "Something went wrong in the construction of FIAT element from UFL element." + \
    #           "Shapes are %s and %s." % (element.value_shape(), ufl_element.value_shape()))

    return element

def create_quadrature(shape, num_points):
    """
    Generate quadrature rule (points, weights) for given shape with
    num_points points in each direction.
    """

    if isinstance(shape, int) and shape == 0:
        return ([()], array([1.0,]))

    if shape in cellname2dim and cellname2dim[shape] == 0:
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

    # Extract the geometric dimension of the points we want to map
    dim = len(points[0]) + 1

    # Special case, don't need to map coordinates on vertices
    if dim == 1:
        return [[(0.0,), (1.0,)][facet]]

    # Get the FIAT reference cell for this dimension
    fiat_cell = reference_cell(dim)

    # Extract vertex coordinates from cell and map of facet index to
    # indicent vertex indices
    vertex_coordinates = fiat_cell.get_vertices()
    facet_vertices = fiat_cell.get_topology()[dim-1]

    #vertex_coordinates = \
    #    {1: ((0.,), (1.,)),
    #     2: ((0., 0.), (1., 0.), (0., 1.)),
    #     3: ((0., 0., 0.), (1., 0., 0.),(0., 1., 0.), (0., 0., 1))}

    # Facet vertices
    #facet_vertices = \
    #    {2: ((1, 2), (0, 2), (0, 1)),
    #     3: ((1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2))}

    # Compute coordinates and map the points
    coordinates = [vertex_coordinates[v] for v in facet_vertices[facet]]
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
        restriction = ufl_element.cell_restriction()
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
    domain = ufl_element.cell_restriction()

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
    #elif isinstance(domain, ufl.Measure):
    #    indices = []
    #    entity_dofs = element.entity_dofs()
    #    for dim, entities in entity_dofs.items():
    #        for entity, index in entities.items():
    #            indices += index
    #    return indices

    else:
        error("Restriction to domain: %s, is not supported." % repr(domain))

