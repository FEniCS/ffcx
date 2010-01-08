__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-11-27"
__copyright__ = "Copyright (C) 2007-2010 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian B. Oelgaard, 2009
# Last changed: 2010-01-08

# Python modules.
from numpy import array

# FIXME: KBO: Use modules from NEW
# FIAT modules.
from FIAT.quadrature import make_quadrature as fiat_make_quadrature
from FIAT.shapes import LINE
from FIAT.shapes import TRIANGLE
from FIAT.shapes import TETRAHEDRON

# FIXME: KBO: Move to fiatinterface
# Glue dictionary
ufl2fiat_shape = {"interval":LINE, "triangle":TRIANGLE, "tetrahedron":TETRAHEDRON}

# FFC modules.
from log import error

def make_quadrature(shape, n, quad_rule=None):
    """
    Generate quadrature rule (points, weights) for given shape with n
    points in each direction. The quadrature rule is generated using
    FIAT and then transformed and scaled to the UFC element.
    """

    # FIXME: Quadrature for vertices (shape == None)
    if shape is None or shape is "vertex":
        return ([()], array([1.0,]))

    # Set scaling and transform
    if shape == "interval":
        offset = array((1.0,))
        scaling = 0.5
    elif shape == "triangle":
        offset = array((1.0, 1.0))
        scaling = 0.25
    elif shape == "tetrahedron":
        offset = array((1.0, 1.0, 1.0))
        scaling = 0.125
    else:
        error("Unknown shape")

    # Get quadrature from FIAT
    q = fiat_make_quadrature(ufl2fiat_shape[shape], n)
    points = q.get_points()
    weights = q.get_weights()

    # Scale from FIAT reference cell to UFC reference cell
    for i in range(len(points)):
        points[i] = tuple(0.5*(array(points[i]) + offset))
        weights[i] *= scaling

    return (points, weights)

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
