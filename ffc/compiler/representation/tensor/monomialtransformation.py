"Transformation of monomial representations of UFL forms."

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2009-03-06 -- 2009-03-06"
__copyright__ = "Copyright (C) 2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# UFL modules
from ufl import BasisFunction, Function

# FFC common modules
from ffc.common.log import ffc_assert

# FFC tensor representation modules
from monomialextraction import MonomialForm

class MonomialIndex:
    pass

class MonomialTransform:
    pass

class MonomialDeterminant:
    pass
    
def transform_monomial_form(monomial_form):
    "Transform monomial form to reference element."

    # Check that we get a Form
    ffc_assert(isinstance(monomial_form, MonomialForm), "Expecting a MonomialForm.")

    # Transform each monomial
    for (integrand, measure) in monomial_form:
        for monomial in integrand.monomials:

            # Reset data for monomial
            determinant = MonomialDeterminant()
            coefficients = []            
            transforms = []
            basis_functions = []

            # Iterate over factors
            for v in monomial.factors:

                # Extract coefficients
                if isinstance(v, Function):
                    coefficient = MonomialCoefficient()
                    coefficients.append(coefficient)

                # Extract transforms
                for d in v.derivative:
                    transform = MonomialTransform()
                    transforms.append(transform)

                # Extract basis functions
                if isinstance(v, BasisFunction):
                    basis_function = MonomialBasisFunction()
                    basis_functions.append(basis_function)
                elif isinstance(v, Function):
                    basis_function = MonomialBasisFunction()
                    basis_functions.append(basis_function)

            # Attach data to monomial
            monomial.determinant = determinant
            monomial.coefficients = coefficients
            monomial.transforms = transforms
            monomial.basis_functions = basis_functions

            # Remove factors
            monomial.factors = []
    
    print "Done"
