__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-11-26 -- 2008-02-11"
__copyright__ = "Copyright (C) 2007-2008 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Marie E. Rognes 2007
# Modified by Kristian B. Oelgaard 2008

from numpy import array

# FFC common modules
from ffc.common.log import error

LINE = 1
TRIANGLE = 2
TETRAHEDRON = 3

coordinates = {LINE: ((0.,), (1.,)),
               TRIANGLE: ((0., 0.), (1., 0.), (0., 1.)),
               TETRAHEDRON: ((0., 0., 0.), (1., 0., 0.),(0., 1., 0.),
                             (0., 0., 1))
               }

facet_identifiers = {TRIANGLE: ((1, 2), (0, 2), (0, 1)),
                     TETRAHEDRON: ((1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2)),
                     }


def get_vertex_coordinates(dim):
    return coordinates[dim]

def get_facet_vertices(dim):
    return facet_identifiers[dim]

def map_to_facet(points, facet):
    """ Mapping of points from the (UFC) reference simplex of
    dimension d - 1 to a given facet on the (UFC) reference simplex of
    dimension d"""
    dim = len(points[0]) + 1

    if dim in [2, 3]:
        facet_vertices = get_facet_vertices(dim)
        vertex_coordinates = get_vertex_coordinates(dim)
    elif dim == 1:
        # Don't need to map coordinates on vertices
        vertex_coordinates = [(0.0,), (1.0,)]
        return [vertex_coordinates[facet]]
    else:
        error("Unable to map points to facet for shape of dimension %d" % dim)

    coordinates = [vertex_coordinates[v] for v in facet_vertices[facet]]

    new_points = []
    for point in points:
        w = (1.0 - sum(point),) + point
        x = tuple(sum([w[i]*array(coordinates[i]) for i in range(len(w))]))
        new_points += [x]
        
    return new_points
