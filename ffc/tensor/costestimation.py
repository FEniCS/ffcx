__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2010-01-25"
__copyright__ = "Copyright (C) 2010 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# FFC modules
from ffc.log import debug

# FFC tensor representation modules
from ffc.tensor.monomialextraction import extract_monomial_form
from ffc.tensor.monomialtransformation import transform_monomial_form

def estimate_cost(integral):
    """
    Estimate cost of tensor representation for integrand. The cost is
    computed as the sum of the number of coefficients and derivatives,
    if the integrand can be represented as a monomial, and -1 if not.
    """

    # Extract monomial integrand
    try:
        monomial_form = extract_monomial_form([integral])
        transform_monomial_form(monomial_form)
    except Exception, exception:
        debug("Monomial extraction failed: " + str(exception))
        return -1

    # Check that we get just one integral
    if not len(monomial_form.integrals) == 1:
        error("Expecting just one integrand.")

    # Compute cost
    (integrand, measure) = monomial_form.integrals[0]
    cost = 0
    for monomial in integrand.monomials:
        cost = max(cost, len(monomial.coefficients) + len(monomial.transforms))

    return cost
