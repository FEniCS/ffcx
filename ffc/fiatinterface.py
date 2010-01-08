__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl) and Anders Logg (logg@simula.no)"
__date__ = "2009-03-06"
__copyright__ = "Copyright (C) 2009 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Garth N. Wells, 2009.
# Modified by Marie Rognes, 2009-2010.
# Last changed: 2010-01-08

# UFL modules
from ufl import FiniteElement as UFLFiniteElement
from ufl import MixedElement as UFLMixedElement
from ufl import ElementRestriction as UFLElementRestriction
from ufl import TensorElement as UFLTensorElement

# FIAT modules
#from FIAT.shapes import LINE, TRIANGLE, TETRAHEDRON
from FIAT_NEW.reference_element import ufc_simplex
from FIAT_NEW.lagrange import Lagrange
from FIAT_NEW.brezzi_douglas_marini import BrezziDouglasMarini
from FIAT_NEW.discontinuous_lagrange import DiscontinuousLagrange
from FIAT_NEW.quadrature import make_quadrature

# FFC modules
from log import debug, error

# FFC fem modules
from quadratureelement import QuadratureElement as FFCQuadratureElement

# Cache for computed elements
_cache = {}

# Mapping from domain to dimension
domain2dim = {"vertex": 0,
              "interval": 1,
              "triangle": 2,
              "tetrahedron": 3}

# Mapping from family name to class
family2class = {"Lagrange": Lagrange,
                "Brezzi-Douglas-Marini": BrezziDouglasMarini,
                "Discontinuous Lagrange": DiscontinuousLagrange}

# Mapping from dimension to number of mesh sub-entities:
entities_per_dim = {1: [2, 1], 2: [3, 3, 1], 3: [4, 6, 4, 1]}

def reference_cell(dim):
    if isinstance(dim, int):
        return ufc_simplex(dim)
    else:
        return ufc_simplex(domain2dim[dim])

def create_element(ufl_element):

    # Use element from cache if in cache
    if ufl_element in _cache:
        debug("Reusing element from cache")
        return _cache[ufl_element]

    # FIXME: hack to avoid circular importing
    from mixedelement import MixedElement as FFCMixedElement

    # Create element
    if isinstance(ufl_element, UFLMixedElement):
        element = FFCMixedElement(ufl_element)
    else:
        element = create_fiat_element(ufl_element)

    # Store in cache
    _cache[ufl_element] = element

    return element

def create_fiat_element(ufl_element):
    "Create FIAT element corresponding to given UFL finite element."

    ElementClass = family2class[ufl_element.family()]
    cell = reference_cell(ufl_element.cell().domain())
    element = ElementClass(cell, ufl_element.degree())

    return element

def create_quadrature(shape, num_points):
    """
    Generate quadrature rule (points, weights) for given shape with
    num_points points in each direction.
    """
    quad_rule = make_quadrature(reference_cell(shape), num_points)
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
    if dim == 1: return [(0.0,), (1.0,)][facet]

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
        w = (1.0 - sum(point),) + point
        x = tuple(sum([w[i]*array(coordinates[i]) for i in range(len(w))]))
        new_points += [x]

    return new_points
