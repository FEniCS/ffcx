__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-10-04"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# FIAT modules
from FIAT.base import quadrature
from FIAT.base import shapes

# FFC modules
from algebra import *
from integrand import *

class Integrator:
    
    """This class is responsible for integrating products of basis
    functions or derivatives of basis functions. The actual work is
    done by FIAT and this class serves as the interface to FIAT
    quadrature rules from within FFC.

    When an Integrator is created, an appropriate quadrature rule is
    chosen based on the total degree of the given product."""

    def __init__(self, basisfunctions):
        "Create Integrator."

        # Check that all functions are defined on the same shape
        for i in range(len(basisfunctions) - 1):
            s0 = basisfunctions[i].element.fiat_shape
            s1 = basisfunctions[i + 1].element.fiat_shape
            if not s0 == s1:
                raise RuntimeError, "BasisFunctions defined on different shapes."

        # All shapes are the same, so pick the first one
        self.fiat_shape = basisfunctions[0].element.fiat_shape

        # Determine the number of required quadrature points
        # FIXME: This should be based on the total degree
        m = 3

        # Create quadrature rule
        self.fiat_quadrature = quadrature.make_quadrature(self.fiat_shape, m)

        # Compensate for different choice of reference cells in FIAT
        # FIXME: Convince Rob and Matt to change their reference cells :-)
        if self.fiat_shape == shapes.TRIANGLE:
            self.vscaling = 0.25  # Area 1/2 instead of 2
            self.dscaling = 2.0   # Scaling of derivative
        elif self.fiat_shape == shapes.TETRAHEDRON:
            self.vscaling = 0.125 # Volume 1/6 instead of 4/3
            self.dscaling = 2.0   # Scaling of derivative
        return

    def __call__(self, basisfunctions, iindices, aindices, bindices):
        "Evaluate integral of product."
        p = Integrand(basisfunctions, iindices, aindices, bindices, self.vscaling, self.dscaling)
        if p.iszero():
            return 0.0 # Makes it a little faster
        else:
            return self.fiat_quadrature(p)
