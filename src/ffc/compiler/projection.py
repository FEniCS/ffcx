__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2005-11-07 -- 2005-11-07"
__copyright__ = "Copyright (c) 2005 Anders Logg"
__license__  = "GNU GPL Version 2"

# Python modules
import sys
import Numeric
import LinearAlgebra

# FIAT modules
from FIAT.quadrature import *

# FFC common modules
sys.path.append("../../")
from ffc.common.exceptions import *

# FFC compiler modules
#from finiteelement import *
#from algebra import *

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
            raise FormError, (function), "Cannot compute projection of given expression."

        # Compute the projection matrix
        P = self.__compute_projection(function.element)

    def __compute_projection(self, e0):
        "Compute projection matrix from e0 to e1."

        # Check if we already know the projection
        name = e0.__repr__()
        if name in self.projections:
            return self.projections[name]

        # Check that the two elements are defined on the same shape
        e1 = self.element
        if not e0.fiat_shape == e1.fiat_shape:
            raise FormError, ((e0, e1)), "Incompatible finite elements for projection."

        # Check that the two elements have the same rank (and rank must be 0 or 1)
        if e0.rank() > 0 or e1.rank() > 0:
            if not (e0.rank() == 1 and e1.rank() == 1):
                raise FormError, ((e0, e1)), "Incompatible element ranks for projection."
            if not (e0.tensordim(0) == e1.tensordim(0)):
                raise FormError, ((e0, e1)), "Incompatible element ranks for projection."

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
        m = e1.spacedim()
        n = e0.spacedim()

        # Create zero order tuple for tables
        dindex = tuple(Numeric.zeros(e0.shapedim()))

        # Compute matrix Q = (vi, vj) for vi in V1
        Q = Numeric.zeros((m, m), Numeric.Float)
        if e0.rank() == 0:
            for i in range(m):
                vi = t1[0][dindex][i]
                for j in range(m):
                    vj = t1[0][dindex][j]
                    sum = 1.0
                    for k in range(len(points)):
                        sum = weights[k] * vi[k] * vj[k]
                    Q[i][j] = sum

        # Compute matrix P = (v_i, w_j) for v_i in V_1 and w_j in V_0
        P = Numeric.zeros((m, n), Numeric.Float)
        for i in range(m):
            for j in range(n):
                P[i][j] = 0.0

        # Compute projection matrix Q^-1 P
        print LinearAlgebra.inverse(Q)
        P = []

        return P

if __name__ == "__main__":

    from finiteelement import *
    from algebra import *
    
    P0 = FiniteElement("Discontinuous Lagrange", "tetrahedron", 0)
    P1 = FiniteElement("Lagrange", "tetrahedron", 1)

    f0 = Function(P0)
    f1 = Function(P1)
    
    pi0 = Projection(P0)
    pi1 = Projection(P1)

#    print f0
#    print f1
    
#    print pi0(f0)
#    print pi0(f1)

#    print pi1(f0)
#    print pi1(f1)

    pi1(f1)
