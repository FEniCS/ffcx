# Copyright (C) 2009-2017 Kristian B. Oelgaard and Martin Sandve Alnæs
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import collections

import numpy

import ufl
from ffcx.ir.representationutils import create_quadrature_points_and_weights
from ufl.classes import Integral
from ufl.sorting import sorted_expr_sum


def compute_quadrature_rules(rules, integral_type, cell):
    """Compute points and weights for a set of quadrature rules."""
    quadrature_rules = {}
    quadrature_rule_sizes = {}
    for scheme, degree in rules:
        # Compute quadrature points and weights
        (points, weights) = create_quadrature_points_and_weights(integral_type, cell, degree,
                                                                 scheme)

        if points is not None:
            points = numpy.asarray(points)

        if weights is None:
            # For custom integrals, there are no points
            num_points = None
        else:
            num_points = len(weights)

        # Assuming all rules with the same number of points are the same
        if num_points in quadrature_rules:
            assert quadrature_rules[num_points][0] == points
            assert quadrature_rules[num_points][0] == weights
            raise RuntimeError(
                "This number of points is already present in the weight table:\n  {}".format(
                    quadrature_rules))

        quadrature_rules[num_points] = (points, weights)
        quadrature_rule_sizes[(scheme, degree)] = num_points

    return quadrature_rules, quadrature_rule_sizes


def accumulate_integrals(itg_data, quadrature_rule_sizes):
    """Group and accumulate integrals according to the number of quadrature points in their rules."""
    if not itg_data.integrals:
        return {}

    # Group integrands by quadrature rule
    sorted_integrands = collections.defaultdict(list)
    if itg_data.integral_type in ufl.custom_integral_types:
        # Should only be one size here, ignoring irrelevant metadata and parameters
        num_points, = quadrature_rule_sizes.values()
        for integral in itg_data.integrals:
            sorted_integrands[num_points].append(integral.integrand())
    else:
        default_scheme = itg_data.metadata["quadrature_degree"]
        default_degree = itg_data.metadata["quadrature_rule"]
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
            sorted_expr_sum(integrands), itg_data.integral_type, itg_data.domain,
            itg_data.subdomain_id, {}, None)
        for num_points, integrands in list(sorted_integrands.items())
    }

    return sorted_integrals
