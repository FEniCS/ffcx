__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-10-04"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# FIAT modules
from FIAT import quadrature
from FIAT import shapes

# FFC common modules
from ffc.common.debug import *

# FFC compiler modules
from algebra import *
from integrand import *

def degree(basisfunctions):
    """Compute total degree of the product of basis functions and
    derivatives of basis functions."""
    q = 0
    for v in basisfunctions:
        q += v.element.degree
        for d in v.derivatives:
            q -= 1
    return q
    
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

        # Determine the number of required quadrature points based on the total
        # degree of the product of basis functions and derivatives of basis functions
        q = degree(basisfunctions)
        m = (q + 1 + 1) / 2 # integer division gives 2m - 1 >= q
        debug("Total degree is %d, using %d quadrature point(s) in each dimension" % (q, m), 1)

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
        v = Integrand(basisfunctions, iindices, aindices, bindices, self.vscaling, self.dscaling)
        if v.iszero():
            debug("      integral is zero, using shortcut", 2)
            return 0.0 # Makes it a little faster
        else:
            return self.fiat_quadrature(v)
