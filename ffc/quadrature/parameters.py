# -*- coding: utf-8 -*-
"Quadrature representation class for UFL"

# Copyright (C) 2009-2014 Kristian B. Oelgaard
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
# Modified by Anders Logg 2009, 2014
# Modified by Martin Sandve Alnæs 2013-2017

# UFL modules
from ufl import custom_integral_types

# FFC modules
from ffc.log import warning


def default_optimize_parameters():
    return {
        "eliminate zeros": False,
        "optimisation": False,
        "ignore ones": False,
        "remove zero terms": False,
        "ignore zero tables": False,
        }


def parse_optimise_parameters(parameters, itg_data):

    # Initialize parameters
    optimise_parameters = default_optimize_parameters()

    # Set optimized parameters
    if parameters["optimize"] and itg_data.integral_type in custom_integral_types:
        warning("Optimization not available for custom integrals, skipping optimization.")

    elif parameters["optimize"]:
        optimise_parameters["ignore ones"] = True
        optimise_parameters["remove zero terms"] = True
        optimise_parameters["ignore zero tables"] = True

        # Do not include this in below if/else clause since we want to
        # be able to switch on this optimisation in addition to the
        # other optimisations.
        if "eliminate_zeros" in parameters:
            optimise_parameters["eliminate zeros"] = True

        if "simplify_expressions" in parameters:
            optimise_parameters["optimisation"] = "simplify_expressions"
        elif "precompute_ip_const" in parameters:
            optimise_parameters["optimisation"] = "precompute_ip_const"
        elif "precompute_basis_const" in parameters:
            optimise_parameters["optimisation"] = "precompute_basis_const"
        # The current default optimisation (for -O) is equal to
        # '-feliminate_zeros -fsimplify_expressions'.
        else:
            # If '-O -feliminate_zeros' was given on the command line,
            # do not simplify expressions
            if "eliminate_zeros" not in parameters:
                optimise_parameters["eliminate zeros"] = True
                optimise_parameters["optimisation"] = "simplify_expressions"

    return optimise_parameters
