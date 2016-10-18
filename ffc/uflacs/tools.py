# -*- coding: utf-8 -*-

# Copyright (C) 2009-2016 Kristian B. Oelgaard and Martin Sandve Alnæs
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

import numbers
import collections
import numpy

from ufl.sorting import sorted_expr_sum
from ufl import custom_integral_types
from ufl.classes import Integral
from ffc.representationutils import create_quadrature_points_and_weights


def collect_quadrature_rules(integrals, default_scheme, default_degree):
    "Collect quadrature rules found in list of integrals."
    rules = set()
    for integral in integrals:
        md = integral.metadata() or {}
        scheme = md.get("quadrature_rule", default_scheme)
        degree = md.get("quadrature_degree", default_degree)
        rule = (scheme, degree)
        rules.add(rule)
    return rules


def compute_quadrature_rules(itg_data):
    "Compute points and weights for a set of quadrature rules."
    # Collect which quadrature rules occur in integrals
    default_scheme = itg_data.metadata["quadrature_degree"]
    default_degree = itg_data.metadata["quadrature_rule"]
    rules = collect_quadrature_rules(
        itg_data.integrals, default_scheme, default_degree)

    quadrature_rules = {}
    quadrature_rule_sizes = {}
    for rule in rules:
        scheme, degree = rule

        # Compute quadrature points and weights
        (points, weights) = create_quadrature_points_and_weights(
            itg_data.integral_type, itg_data.domain.ufl_cell(), degree, scheme)

        if points is not None:
            points = numpy.asarray(points)

        if weights is None:
            # For custom integrals, there are no points
            num_points = None
        else:
            num_points = len(weights)

        # Assuming all rules with the same number of points are equal
        if num_points in quadrature_rules:
            assert quadrature_rules[num_points][0] == points
            assert quadrature_rules[num_points][0] == weights
            error("This number of points is already present in the weight table:\n  %s" % (quadrature_rules,))

        quadrature_rules[num_points] = (points, weights)
        quadrature_rule_sizes[rule] = num_points
    return quadrature_rules, quadrature_rule_sizes


def accumulate_integrals(itg_data, quadrature_rule_sizes):
    """Group and accumulate integrals according to the number
    of quadrature points in their rules.
    """
    if not itg_data.integrals:
        return {}

    default_scheme = itg_data.metadata["quadrature_degree"]
    default_degree = itg_data.metadata["quadrature_rule"]

    # Group integrands by quadrature rule
    sorted_integrands = collections.defaultdict(list)
    for integral in itg_data.integrals:
        md = integral.metadata() or {}
        scheme = md.get("quadrature_rule", default_scheme)
        degree = md.get("quadrature_degree", default_degree)
        rule = (scheme, degree)
        num_points = quadrature_rule_sizes[rule]
        sorted_integrands[num_points].append(integral.integrand())

    # Accumulate integrands in a canonical ordering defined by UFL
    sorted_integrals = {
        num_points: Integral(
            sorted_expr_sum(integrands),
            itg_data.integral_type,
            itg_data.domain,
            itg_data.subdomain_id,
            {},
            None)
        for num_points, integrands in list(sorted_integrands.items())
        }
    return sorted_integrals
