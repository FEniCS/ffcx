__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-10-04"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# FIAT modules
from FIAT.base import quadrature

# FFC modules
from algebra import *

class Integrator:
    
    """This class is responsible for integrating products of basis
    functions or derivatives of basis functions. The actual work is
    done by FIAT and this class serves as the interface to FIAT
    quadrature rules from within FFC.

    When an Integrator is created, an appropriate quadrature rule is
    chosen based on the total degree of the given Product."""

    def __init__(self, product):
        "Create Integrator."

        # Check that all functions are defined on the same shape
        for i in range(len(product.factors) - 1):
            s0 = product.factors[i].basisfunction.element.fiat_shape
            s1 = product.factors[i + 1].basisfunction.element.fiat_shape
            if not s0 == s1:
                raise RuntimeError, "BasisFunctions defined on different shapes."

        # All shapes are the same, so pick the first one
        self.fiat_shape = product.factors[0].basisfunction.element.fiat_shape

        # Determine the number of required quadrature points
        # FIXME: This should be based on the total degree
        m = 4

        # Create quadrature rule
        self.fiat_quadrature = quadrature.make_quadrature(self.fiat_shape, m)

        return

    def __call__(self, product, indices0, indices1):
        "Evaluate integral of Product."
        p = ProductFunction(product, indices0, indices1)
        return self.fiat_quadrature(p)

class ProductFunction:
    
    """This class wraps around a Product to create a callable
    object. We put this functionality here rather than in the Product
    class itself to avoid putting to many things into the algebra
    module."""

    def __init__(self, product, indices0, indices1):
        "Create ProductFunction."
        if not isinstance(product, Product):
            raise RuntimeError, "Product expected."
        self.product = product
        self.indices0 = indices0
        self.indices1 = indices1
        return

    def __call__(self, x):
        "Evaluate Product at given point."
        tmp = 1.0
        for f in self.product.factors:
            tmp = tmp * self.__eval(f, x)
        return tmp

    def __eval(self, factor, x):
        "Evaluate Factor at given point."

        # Get basis function
        i = factor.basisfunction.index(self.indices0, self.indices1, [])
        v = factor.basisfunction.element.basis[i]

        # Differentiate the basis function
        for d in factor.derivatives:
            i = d.index(self.indices0, self.indices1, [])
            v = v.deriv(i)

        # Evaluate basis function
        return v(x)
