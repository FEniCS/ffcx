__author__ = "Garth N. Wells (gnw20@cam.ac.uk)"
__date__ = "2011-04-19"
__copyright__ = "Copyright (C) 2011 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# FIAT modules
import ufl
import FIAT

# FFC modules
from ffc.log import debug, error
from ffc.fiatinterface import reference_cell

# Dictionary mapping from domain (cell) to dimension
from ufl.geometry import domain2dim

def create_quadrature(shape, num_points):
    """
    Generate quadrature rule (points, weights) for given shape with
    num_points points in each direction.
    """

    # FIXME: KBO: Can this be handled more elegantly?
    # Handle point case
    if isinstance(shape, int) and shape == 0 or domain2dim[shape] == 0:
        return ([()], array([1.0,]))

    quad_rule = FIAT.make_quadrature(reference_cell(shape), num_points)
    return quad_rule.get_points(), quad_rule.get_weights()

