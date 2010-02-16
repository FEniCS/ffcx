__author__ = "Kristian B. Oelgaard (k.b.oelgaard@gmail.com) and Anders Logg (logg@simula.no)"
__date__ = "2009-03-06"
__copyright__ = "Copyright (C) 2009-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Garth N. Wells, 2009.
# Modified by Marie Rognes, 2009-2010.
# Last changed: 2010-02-01

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

# Cache for computed elements
_cache = {}

# FIXME: KBO: Should stuff like, domain2dim and entities_per_dim be in UFC
# instead? The same goes for similar dictionaries in UFL (geometry.py). After
# all both FFC and UFL complies with UFC or not?
# Mapping from domain to dimension
# FIXME: AL: They should be in UFL (and probably are there already). They
# FIXME: can't be in UFC since UFL cannot depend on UFC.
domain2dim = {"vertex": 0,
              "interval": 1,
              "triangle": 2,
              "tetrahedron": 3}

# Mapping from dimension to number of mesh sub-entities:
entities_per_dim = {1: [2, 1], 2: [3, 3, 1], 3: [4, 6, 4, 1]}

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

    # Create restricted element(implemented by FFC)
    elif isinstance(ufl_element, ufl.ElementRestriction):
        element = _create_restricted_element(ufl_element)

    else:
        error("Cannot handle this element type: %s" % str(ufl_element))

    # Store in cache
    _cache[ufl_element] = element

    return element

def _create_fiat_element(ufl_element):
    "Create FIAT element corresponding to given finite element."

    family = ufl_element.family()

    # FIXME: AL: Should this really be here?
    # Handle QuadratureElement
    if family == "Quadrature":
        return FFCQuadratureElement(ufl_element)

    # Check if finite element family is supported by FIAT
    elif not family in FIAT.element_classes:
        error("Sorry, finite element of type \"%s\" are not supported by FIAT.", family)

    # Create FIAT finite element
    ElementClass = FIAT.element_classes[family]
    cell = reference_cell(ufl_element.cell().domain())
    element = ElementClass(cell, ufl_element.degree())

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

    if domain:
        ufl_element = ufl.ElementRestriction(ufl_element, domain)

    elements += [create_element(ufl_element)]
    return elements

def _create_restricted_element(ufl_element):
    "Create an FFC representation for an UFL ElementRestriction."

    if not isinstance(ufl_element, ufl.ElementRestriction):
        error("create_restricted_element expects an ufl.ElementRestriction")

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

def _indices(element, domain):
    "Extract basis functions indices that correspond to domain."

    if isinstance(domain, ufl.Cell):
        dim = domain.topological_dimension()
        entity_dofs = element.entity_dofs()
        indices = []
        # FIXME: KBO: This will only make it possible to restrict up to a given
        # topological dimension. What if one wants to restrict to the dofs on
        # the interior of a cell?
        for dim in range(domain.topological_dimension() + 1):
            entities = entity_dofs[dim]
            for (entity, index) in entities.iteritems():
                indices += index
        return indices
    # Just extract all indices to make handling in RestrictedElement uniform.
    elif isinstance(domain, ufl.Measure):
        indices = []
        entity_dofs = element.entity_dofs()
        for dim, entities in entity_dofs.items():
            for entity, index in entities.items():
                indices += index
        return indices
    else:
        error("Restriction to domain: %s, is not supported." % repr(domain))

