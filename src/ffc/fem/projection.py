__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2005-11-07 -- 2007-01-23"
__copyright__ = "Copyright (C) 2005-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# Modified by Garth N. Wells 2006

# Python modules
import sys
import numpy
import numpy.linalg

# FIAT modules
from FIAT.quadrature import *

# FFC common modules
sys.path.append("../../")
from ffc.common.exceptions import *

# FFC compiler.language modules
from ffc.compiler.language.algebra import *

class Projection:
    """A Projection represents the local L2 projection onto a given
    finite element.

    Attributes:

        element     - the FiniteElement defining the projection
        projections - dictionary of cached projection matrices
        
    """

    def __init__(self, element):
        "Create Projection onto given FiniteElement."

        self.element = element
        self.projections = {}

    def __call__(self, function):
        "Compute projection of given Function."

        # Check that we got a Function
        if not isinstance(function, Function):
            raise FormError, (function, "Projections are only supported for Functions.")

        # Check that we have not already computed the projection
        if not function.P == None:
            raise FormError, (function, "Only one projection can be applied to each Function.")

        # Compute the projection matrix
        P = self.__compute_projection(function.e0)

        # Create new function and set projection
        f = Function(function)
        f.n1 = Index("projection")
        f.e1 = self.element
        f.P  = P

        return f

    def __compute_projection(self, e0):
        "Compute projection matrix from e0 to e1."

        # Check if we already know the projection
        name = e0.signature()
        if name in self.projections:
            print "Projection cached, reusing projection"
            return self.projections[name]

        # Check that the two elements are defined on the same shape
        e1 = self.element
        if not e0.fiat_shape == e1.fiat_shape:
            raise FormError, (((e0, e1)), "Incompatible finite elements for projection.")

        # Check that the two elements have the same rank (and rank must be 0 or 1)
        if e0.value_rank() > 0 or e1.value_rank() > 0:
            if not (e0.value_rank() == 1 and e1.value_rank() == 1):
                raise FormError, (((e0, e1)), "Incompatible element ranks for projection.")
            if not (e0.value_dimension(0) == e1.value_dimension(0)):
                raise FormError, (((e0, e1)), "Incompatible element ranks for projection.")
            vectordim = e0.value_dimension(0) # Both dimensions are equal
        rank = e0.value_rank() # Both ranks are equal

        # Get quadrature points and weights for integrals
        q = e1.degree() + max(e0.degree(), e1.degree()) # same points for both integrals
        m = (q + 1 + 1) / 2 # integer division gives 2m - 1 >= q
        quadrature = make_quadrature(e0.fiat_shape, m)
        points = quadrature.get_points()
        weights = quadrature.get_weights()

        # Tabulate basis functions at quadrature points
        t0 = e0.tabulate(0, points)
        t1 = e1.tabulate(0, points)

        # Get space dimensions of V0 and V1
        m = e1.space_dimension()
        n = e0.space_dimension()

        # Create zero order tuple for tables
        dindex = tuple(numpy.zeros(e0.shape_dimension(), dtype = numpy.int))

        # Compute matrix Q = (vi, vj) for vi in V1
        Q = numpy.zeros((m, m), dtype = numpy.float)
        if rank == 0:
            for i in range(m):
                vi = t1[0][dindex][i]
                for j in range(m):
                    vj = t1[0][dindex][j]
                    sum = 0.0
                    for k in range(len(points)):
                        sum += weights[k] * vi[k] * vj[k]
                    Q[i][j] = sum
        else:
            for l in range(vectordim):
                for i in range(m):
                    vi = t1[l][0][dindex][i]
                    for j in range(m):
                        vj = t1[l][0][dindex][j]
                        sum = 0.0
                        for k in range(len(points)):
                            sum += weights[k] * vi[k] * vj[k]
                        Q[i][j] += sum

        # Compute matrix P = (v_i, w_j) for v_i in V_1 and w_j in V_0
        P = numpy.zeros((m, n), dtype = numpy.float)
        if rank == 0:
            for i in range(m):
                vi = t1[0][dindex][i]
                for j in range(n):
                    wj = t0[0][dindex][j]
                    sum = 0.0
                    for k in range(len(points)):
                        sum += weights[k] * vi[k] * wj[k]
                    P[i][j] = sum
        else:
            for l in range(vectordim):
                for i in range(m):
                    vi = t1[l][0][dindex][i]
                    for j in range(n):
                        wj = t0[l][0][dindex][j]
                        sum = 0.0
                        for k in range(len(points)):
                            sum += weights[k] * vi[k] * wj[k]
                        P[i][j] += sum

        # Compute projection matrix Q^-1 P
        P = numpy.dot(numpy.linalg.inv(Q), P)

        # Save projection matrix for later so it can be reused
        self.projections[name] = P
        
        return P
