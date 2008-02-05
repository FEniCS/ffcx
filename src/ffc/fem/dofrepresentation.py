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
    """A degree of freedom is represented by its points, direction
    and weights.

    name            = An explanatory name for the functional
    points          = Points for which the functional is to be evaluated at.
    direction       = Direction for functional components 
    weights         = Values used for integral moments.

    Given a function f, the DofRepresentation represents the value:

       return weights[i]*f[k](points[i])*direction[k] 

    """

    def __init__(self, dof, points=None, direction=None, weights=None):
        "Create DofRepresentation"

        # Create dof representation from other dof representation:
        # Copy contents.
        if isinstance(dof, DofRepresentation):
            self.name = dof.name
            self.points = [p for p in dof.points]
            self.direction = dof.direction
            self.weights = [w for w in dof.weights]
        else:
            self.name = dof
            if not points:
                print "Raise Error: No points"
            else:
                self.points = points
            if len(direction) < 1:
                self.direction = [1]
            else:
                self.direction = direction
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
        return len(self.weights)

    def value_dim(self):
        return len(self.direction)

    def shift_direction(self, n, shift):
        """Take the direction, add k zeros at the beginning and ensure
        that the new direction has length n"""
        newdirection = [0]*n
        newdirection[shift:shift+self.value_dim()] = self.direction
        self.direction = newdirection

    def pad_points_and_weights(self, n):
        """ For easy c++ representation, pad points and weights with
        zeros at the end. Return the original length of the entities"""
        k = self.cell_dimension()
        m = self.num_of_points()
        ptail = [[0]*k]*(n-m)
        wtail = [0]*(n-m)
        self.points = self.points + ptail
        self.weights = self.weights + wtail
        return m

    def __str__(self):
        """ Pretty print """
        return str(self.name)
