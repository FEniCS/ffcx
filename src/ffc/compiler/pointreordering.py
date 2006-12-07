"""This module implements reordering of quadrature points on
a shape. This is used to reorder the quadrature points on
a facet depending on the alignment of the facet with respect
to a neighboring cell."""

__author__ = "Kristian Oelgaard (k.b.oelgaard@tudelft.nl) and Anders Logg (logg@simula.no)"
__date__ = "2005-09-06 -- 2006-12-06"
__copyright__ = "Copyright (C) 2006 Kristian Oelgaard and Anders Logg"
__license__  = "GNU GPL Version 2"

# FIAT modules
from FIAT.shapes import *

def reorder_points(points, shape, alignment):
    """Reorder points on given reference shape for given alignment."""

    if shape == LINE:
        return __reorder_line(points, alignment)
    elif shape == TRIANGLE:
        return __reorder_triangle(points, alignment)
    else:
        raise RuntimeError, "Unable to reorder quadrature points. Unknown shape: " + str(shape)

def __reorder_line(points, alignment):
    """Reorder points on the (FIAT) reference line."""

    def Phi0(x):
        return 0.5 - 0.5*x[0]

    def Phi1(x):
        return 0.5 + 0.5*x[0]

    if alignment == 0:
        p0 = -1.0
        p1 =  1.0
    elif alignment == 1:
        p0 =  1.0
        p1 = -1.0

    new_points = [(Phi0(x)*p0 + Phi1(x)*p1,) for x in points]

    return new_points
    
def __reorder_triangle(x, alignment):
    """Reorder points on the (FIAT) reference triangle."""
    # FIXME: Write code here
    print "triangle ", x, "align ", alignment
    return x
