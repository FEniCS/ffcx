__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-10-04 -- 2005-09-14"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# FIAT modules
from FIAT.quadrature import *
from FIAT.shapes import *

# FFC common modules
from ffc.common.debug import *
from ffc.common.exceptions import *

# FFC compiler modules
from algebra import *

def degree(basisfunctions):
    """Compute total degree of the product of basis functions and
    derivatives of basis functions."""
    q = 0
    for v in basisfunctions:
        q += v.element.degree()
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
            s0 = basisfunctions[i].element.shape()
            s1 = basisfunctions[i + 1].element.shape()
            if not s0 == s1:
                raise RuntimeError, (basisfunctions, "BasisFunctions defined on different shapes.")

        # All shapes are the same, so pick the first one
        self.shape = basisfunctions[0].element.shape()

        # Initialize quadrature
        self.__init_quadrature(basisfunctions)

        # Tabulate basis functions at quadrature points
        self.__init_table(basisfunctions)
        
        return

    def __init_quadrature(self, basisfunctions):
        "Create quadrature rule."

        # Determine the number of required quadrature points based on the total
        # degree of the product of basis functions and derivatives of basis functions
        q = degree(basisfunctions)
        m = (q + 1 + 1) / 2 # integer division gives 2m - 1 >= q
        debug("Total degree is %d, using %d quadrature point(s) in each dimension" % (q, m), 1)

        # Create quadrature rule and get points and weights
        quadrature = make_quadrature(self.shape, m)
        self.points = quadrature.get_points()
        self.weights = quadrature.get_weights()

        # Compensate for different choice of reference cells in FIAT
        # FIXME: Convince Rob and Matt to change their reference cells :-)
        if self.shape == TRIANGLE:
            self.vscaling = 0.25  # Area 1/2 instead of 2
            self.dscaling = 2.0   # Scaling of derivative
        elif self.shape == TETRAHEDRON:
            self.vscaling = 0.125 # Volume 1/6 instead of 4/3
            self.dscaling = 2.0   # Scaling of derivative

        return
        
    def __init_table(self, basisfunctions):
        "Create table of basisfunctions at quadrature points."

        # Compute maximum number of derivatives for each element
        num_derivatives = {}
        for v in basisfunctions:
            element = v.element
            order = len(v.derivatives)
            if element in num_derivatives:
                if order > num_derivatives[element]:
                    num_derivatives[element] = order
            else:
                num_derivatives[element] = order

        # Call FIAT to tabulate the basis functions for each element
        self.table = {}
        for element in num_derivatives:
            order = num_derivatives[element]
            self.table[element] = element.tabulate(order, self.points)

    def __call__(self, basisfunctions, iindices, aindices, bindices):
        "Evaluate integral of product."

        integrand = self.vscaling * self.weights
        
        for v in basisfunctions:
            (element, vindex, cindex, dorder, dindex) = v(iindices, aindices, bindices)
            if element.rank() > 0:
                value = self.table[element][cindex[0]][dorder][dindex][vindex]
            else:
                value = self.table[element][dorder][dindex][vindex]
            integrand = integrand * value * pow(self.dscaling, dorder)
                
        return sum(integrand)
