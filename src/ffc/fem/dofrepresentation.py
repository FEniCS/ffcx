__author__ = "Marie E. Rognes (meg@cma.uio.no)"
__date__ = "2008-01-22 -- "
__copyright__ = "Copyright (C) 2008 Marie Rognes"
__license__  = "GNU GPL version 3 or any later version"

# Python modules
import sys

# FFC common modules
from ffc.common.debug import *
from ffc.common.utils import *

# FFC fem modules

class DofRepresentation:
    """A degree of freedom is represented by its point(s), direction(s)
    and weight(s).

    name            = An explanatory name for the functional
    points          = Points for which the functional is to be evaluated at.
    directions      = Directions for functional components 
    weights         = Values used for integral moments.

    Given a function f, the DofRepresentation represents the value:

       return weights[i]*f[k](points[i])*directions[i][k] 

    This is intended to cover point values, directional components and
    integral moments.
    """

    def __init__(self, dof, points=None, directions=None, weights=None):
        "Create DofRepresentation"
        # Create dof representation from other dof representation:
        # Copy contents.
        if isinstance(dof, DofRepresentation):
            self.name = dof.name
            self.points = [p for p in dof.points]
            self.directions = [d for d in dof.directions]
            self.weights = [w for w in dof.weights]
        else:
            self.name = dof
            if len(points) < 1:
                points = []
            else:
                self.points = points
            if len(directions) < 1:
                self.directions = [[1]]*self.num_of_points()
            else:
                self.directions = directions
            if len(weights) < 1:
                self.weights = [1]*self.num_of_points()
            else:
                self.weights = weights
        return

    def num_of_points(self):
        return len(self.points)

    def cell_dimension(self):
        "Return the dimension of the point values"
        return pick_first([len(pt) for pt in self.points])

    def num_of_weights(self):
        """Return the number of weights. Should match the number of
        points."""
        return len(self.weights)

    def value_dim(self):
        "Return the value dimension of the directions"
        return pick_first([len(d) for d in self.directions])

    def shift_directions(self, n, shift):
        """Take the directions, add k zeros at the beginning and
        ensure that the new directions has length n"""
        d = self.value_dim()
        for i in range(len(self.directions)):
            newdirection = [0]*n
            newdirection[shift:shift+d] = self.directions[i]
            self.directions[i] = newdirection

    def pad_points_and_weights(self, n):
        """ For easy c++ representation, pad points, directions and
        weights with zeros at the end. Return the original length of
        the entities"""
        k = self.cell_dimension()
        d = self.value_dim()
        m = self.num_of_points()
        ptail = [[0]*k]*(n-m)
        dtail = [[0]*d]*(n-m)
        wtail = [0]*(n-m)
        self.points = self.points + ptail
        self.directions = self.directions + dtail
        self.weights = self.weights + wtail
        return m

    def __str__(self):
        """ Pretty print """
        return str(self.name)
