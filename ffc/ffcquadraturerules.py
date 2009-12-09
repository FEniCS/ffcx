__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-11-27"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian B. Oelgaard, 2009
# Last changed: 2009-12-09

# Python modules.
from numpy import array

# FIAT modules.
from FIAT.quadrature import make_quadrature as fiat_make_quadrature
from FIAT.shapes import LINE
from FIAT.shapes import TRIANGLE
from FIAT.shapes import TETRAHEDRON

# Glue dictionary
ufl2fiat_shape = {"interval":LINE, "triangle":TRIANGLE, "tetrahedron":TETRAHEDRON}

# FFC modules.
from log import error

def make_quadrature(shape, n, quad_rule=None):
    """Generate quadrature rule (points, weights) for given shape
    with n points in each direction. Quadrature rule is generated
    using FIAT and then transformed and scaled to the UFC element"""

    # FIXME: Quadrature for vertices (shape == None)
    if shape is None:
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

if __name__ == "__main__":

    from FIAT.shapes import *

    print make_quadrature(TRIANGLE, 2)
