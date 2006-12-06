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

def reorder_points(x, shape, alignment):
    """Reorder points on given reference shape for given alignment."""

    if shape == LINE:
        return __reorder_triangle(x, alignment)
    elif shape == TRIANGLE:
        return __reorder_tetrahedron(x, alignmnent)
    else:
        raise RuntimeError, "Unable to reorder quadrature points, unknown shape: " + str(shape)

def __reorder_line(x, alignment):
    """Reorder points on the (FIAT) reference triangle."""
    # FIXME: Write code here
    return x

def __reorder_triangle(x, alignment):
    """Reorder points on the (FIAT) reference triangle."""
    # FIXME: Write code here
    return x
