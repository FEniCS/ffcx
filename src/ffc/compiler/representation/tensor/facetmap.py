"""Mapping of points from the (UFC) reference simplex of dimension d - 1
to a given facet on the (UFC) reference simplex of dimension d"""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-11-26 -- 2007-11-26"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

from numpy import array

def map_to_facet(points, facet):
    "Map points to given facet"

    dim = len(points[0]) + 1
    
    if dim == 2:
        facet_vertices = [(1, 2), (0, 2), (0, 1)]
        vertex_coordinates = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
    elif dim == 3:
        facet_vertices = [(1 , 2 , 3), (0, 2, 3), (0, 1, 3), (0, 1, 2)]
        vertex_coordinates = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
    else:
        raise RuntimeError, "Unable to map points to facet for shape of dimension %d" % dim

    coordinates = [vertex_coordinates[v] for v in facet_vertices[facet]]

    new_points = []
    for point in points:
        w = (1.0 - sum(point),) + point
        x = tuple(sum([w[i]*array(coordinates[i]) for i in range(len(w))]))
        new_points += [x]
        
    return new_points
