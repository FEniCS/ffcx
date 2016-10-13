# -*- coding: utf-8 -*-
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
# Last changed: 2014-04-15

# FFC modules
from ffc.log import debug, error

# FFC tensor representation modules
from ffc.tensor.monomialextraction import extract_monomial_form
from ffc.tensor.monomialtransformation import transform_monomial_form

def estimate_cost(integral, function_replace_map):
    """
    Estimate cost of tensor representation for integral. The cost is
    computed as the sum of the number of coefficients and derivatives,
    if the integrand can be represented as a monomial, and -1 if not.
    """

    # Check that integral type is supported
    supported = ["cell", "exterior_facet", "interior_facet"]
    if not integral.integral_type() in supported:
        return -1

    # Extract monomial integrand
    integrand = integral.integrand()
    try:
        monomial_form = extract_monomial_form([integrand], function_replace_map)
        transform_monomial_form(monomial_form)
    except Exception as exception:
        debug("Monomial extraction failed: " + str(exception))
        return -1

    # Check that we get just one integrand
    if not len(monomial_form) == 1:
        error("Expecting just one integrand.")

    # Compute cost
    cost = 0
    for integrand in monomial_form:
        for monomial in integrand.monomials:
            cost = max(cost, len(monomial.coefficients) + len(monomial.transforms))
    return cost
