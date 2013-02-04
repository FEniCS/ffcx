# Copyright (C) 2010 Anders Logg
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.
#
# First added:  2010-01-25
# Last changed: 2010-01-25

# FFC modules
from ffc.log import debug

# FFC tensor representation modules
from ffc.tensor.monomialextraction import extract_monomial_form
from ffc.tensor.monomialtransformation import transform_monomial_form

def estimate_cost(integral, function_replace_map):
    """
    Estimate cost of tensor representation for integrand. The cost is
    computed as the sum of the number of coefficients and derivatives,
    if the integrand can be represented as a monomial, and -1 if not.
    """

    # TODO: Implement this to take integrand instead of integral? May allow simplification inn caller code.

    # Extract monomial integrand
    try:
        monomial_form = extract_monomial_form([integral], function_replace_map)
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
