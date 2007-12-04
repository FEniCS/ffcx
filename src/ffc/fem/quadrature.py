__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-11-27 -- 2007-11-29"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

from numpy import array

from FIAT.quadrature import make_quadrature as fiat_make_quadrature
from FIAT.shapes import LINE, TRIANGLE, TETRAHEDRON

def make_quadrature(shape, n):
    """Generate quadrature rule (points, weights) for given shape
    with n points in each direction. Quadrature rule is generated
    using FIAT and then transformed and scaled to the UFC element"""

    # FIXME: Quadrature for vertices (shape == None)
    if shape:
        # Get quadrature from FIAT
        q = fiat_make_quadrature(shape, n)
        points = q.get_points()
        weights = q.get_weights()

        # FIXME: Temporary until things work
        return (points, weights)
    else:
        return ([()], array([1.0,]))

    # Set scaling and transform
    if shape == LINE:
        offset = array((1.0,))
        scaling = 0.5
    elif shape == TRIANGLE:
        offset = array((1.0, 1.0))
        scaling = 0.25
    elif shape == TETRAHEDRON:
        offset = array((1.0, 1.0, 1.0))
        scaling = 0.125
    else:
        raise RuntimeError, "Unknown shape"
    
    # Scale from FIAT reference cell to UFC reference cell    
    for i in range(len(points)):
        points[i] = tuple(0.5*(array(points[i]) + offset))
        weights[i] *= scaling

    return (points, weights)

if __name__ == "__main__":

    from FIAT.shapes import *

    print make_quadrature(TRIANGLE, 2)
