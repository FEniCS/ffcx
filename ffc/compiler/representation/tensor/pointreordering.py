"""This module implements reordering of quadrature points on facets,
which is necessary when integrating over an interior facet where two
cells meet. The two cells then need to agree on the orientation of
the common facet. This would not be necessary if FIAT were to use
the UFC ordering convention for mesh entities."""

__author__ = "Kristian Oelgaard (k.b.oelgaard@tudelft.nl) and Anders Logg (logg@simula.no)"
__date__ = "2006-12-19 -- 2007-05-14"
__copyright__ = "Copyright (C) 2006-2007 Kristian Oelgaard and Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# FIAT modules
from FIAT.shapes import *

from numpy import array

def reorder_points(points, shape, facet):
    """Reorder points on given reference shape for given facet."""

    if shape == TRIANGLE:
        return __reorder_edge(points, facet)
    elif shape == TETRAHEDRON:
        return __reorder_face(points, facet)
    else:
        raise RuntimeError, "Unable to reorder quadrature points. Unknown shape: " + str(shape)

def __reorder_edge(points, facet):
    """Reorder points on edge of (FIAT) reference triangle."""

    def Phi0(x):
        return 0.5 - 0.5*x[0]

    def Phi1(x):
        return 0.5 + 0.5*x[0]

    if facet == 0:
        p0 = -1.0
        p1 =  1.0
    elif facet == 1:
        # This is the case where FIAT differs from the UFC specification
        p0 =  1.0
        p1 = -1.0
    elif facet == 2:
        p0 = -1.0
        p1 =  1.0
    else:
        raise RuntimeError, "Illegal facet for reordering of points on facet."

    new_points = [(Phi0(x)*p0 + Phi1(x)*p1, ) for x in points]

    return new_points
    
def __reorder_face(points, facet):
    """Reorder points on facet of (FIAT) reference triangle."""
    
    def Phi0(x):
        return -0.5*x[0] - 0.5*x[1]

    def Phi1(x):
        return 0.5*x[0] + 0.5

    def Phi2(x):
        return 0.5*x[1] + 0.5

    if facet == 0:
        p0 = (-1.0, -1.0)
        p1 = (-1.0, 1.0)
        p2 = (1.0, -1.0)
    elif facet == 1:
        p0 = (-1.0, 1.0)
        p1 = (-1.0, -1.0)
        p2 = (1.0, -1.0)
    elif facet == 2:
        p0 = (-1.0, 1.0)
        p1 = (1.0, -1.0)
        p2 = (-1.0, -1.0)
    elif facet == 3:
        # This is the only case where FIAT does NOT differ from the UFC specification
        p0 = (-1.0, -1.0)
        p1 = (1.0, -1.0)
        p2 = (-1.0, 1.0)
    else:
        raise RuntimeError, "Illegal facet for reordering of points on facet."

    p0 = array(p0)
    p1 = array(p1)
    p2 = array(p2)

    new_points = [tuple(Phi0(x)*p0 + Phi1(x)*p1 + Phi2(x)*p2) for x in points]

    return new_points
