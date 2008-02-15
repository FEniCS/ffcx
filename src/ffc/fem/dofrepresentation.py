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
        if isinstance(dof, DofRepresentation):
            self.name = dof.name
            self.points = listcopy(dof.points)
            self.directions = listcopy(dof.directions)
            self.weights = listcopy(dof.weights)
        else:
            self.name = dof
            if points == None:
                # No points indicates that something is not implemented...
                self.points = []
            else:
                self.points = points
            if directions == None:
                self.directions = [(1,)]*len(self.points)
            else:
                self.directions = directions
            if weights == None:
                self.weights = [1]*len(self.points)
            else:
                self.weights = weights

        # Check that the dimensions match:
        assert len(self.points) == len(self.weights) == len(self.directions), \
               "Mismatch points/directions/weights"
        return

    def num_of_points(self):
        return len(self.points)

    def cell_dimension(self):
        "Return the dimension of the point values"
        return pick_first([len(pt) for pt in self.points])

    def num_of_weights(self):
        """Return the number of weights."""
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
        # Add (n-m) number of zero (vectors) at the end
        self.points = self.points + [tuple([0]*k)]*(n-m)
        self.directions = self.directions + [tuple([0]*d)]*(n-m) 
        self.weights = self.weights + [0]*(n-m) 
        return m

    def __str__(self):
        """ Pretty print """
        return str(self.name)
