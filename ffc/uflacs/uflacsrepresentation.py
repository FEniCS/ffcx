# -*- coding: utf-8 -*-
# Copyright (C) 2013-2017 Martin Sandve Alnæs
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import logging

from ffc.fiatinterface import create_element
from ffc.representationutils import initialize_integral_ir
from ffc.uflacs.build_uflacs_ir import build_uflacs_ir
from ffc.uflacs.tools import (accumulate_integrals, collect_quadrature_rules,
                              compute_quadrature_rules)
from ufl import custom_integral_types
from ufl.algorithms import replace
from ufl.utils.sorting import sorted_by_count
import ufl

logger = logging.getLogger(__name__)


def compute_integral_ir(itg_data, form_data, form_id, element_numbers, classnames, parameters):
    """Compute intermediate represention of integral."""

    logger.info("Computing uflacs representation")

    # Initialise representation
    ir = initialize_integral_ir("uflacs", itg_data, form_data, form_id)

    # Store element classnames
    ir["classnames"] = classnames

    # Get element space dimensions
    unique_elements = element_numbers.keys()
    ir["element_dimensions"] = {
        ufl_element: create_element(ufl_element).space_dimension()
        for ufl_element in unique_elements
    }

    # Create dimensions of primary indices, needed to reset the argument 'A'
    # given to tabulate_tensor() by the assembler.
    argument_dimensions = [
        ir["element_dimensions"][ufl_element] for ufl_element in form_data.argument_elements
    ]

    # Compute shape of element tensor
    if ir["integral_type"] == "interior_facet":
        ir["tensor_shape"] = [2 * dim for dim in argument_dimensions]
    else:
        ir["tensor_shape"] = argument_dimensions

    integral_type = itg_data.integral_type
    cell = itg_data.domain.ufl_cell()

    if integral_type in custom_integral_types:
        # Set quadrature degree to twice the highest element degree, to get
        # enough points to identify basis functions via table computations
        max_element_degree = max([1] + [ufl_element.degree() for ufl_element in unique_elements])
        rules = [("default", 2 * max_element_degree)]
        quadrature_integral_type = "cell"
    else:
        # Collect which quadrature rules occur in integrals
        default_scheme = itg_data.metadata["quadrature_degree"]
        default_degree = itg_data.metadata["quadrature_rule"]
        rules = collect_quadrature_rules(itg_data.integrals, default_scheme, default_degree)
        quadrature_integral_type = integral_type

    # Compute actual points and weights
    quadrature_rules, quadrature_rule_sizes = compute_quadrature_rules(
        rules, quadrature_integral_type, cell)

    # Store quadrature rules in format { num_points: (points, weights) }
    ir["quadrature_rules"] = quadrature_rules

    # Store the fake num_points for analysis in custom integrals
    if integral_type in custom_integral_types:
        ir["fake_num_points"], = quadrature_rules.keys()

    # Group and accumulate integrals on the format { num_points: integral data }
    sorted_integrals = accumulate_integrals(itg_data, quadrature_rule_sizes)

    # Build coefficient numbering for UFC interface here, to avoid
    # renumbering in UFL and application of replace mapping

    # Using the mapped coefficients, numbered by UFL
    coefficient_numbering = {}
    sorted_coefficients = sorted_by_count(form_data.function_replace_map.keys())
    for i, f in enumerate(sorted_coefficients):
        g = form_data.function_replace_map[f]
        coefficient_numbering[g] = i
        assert i == g.count()

    # Replace coefficients so they all have proper element and domain for what's to come
    # TODO: We can avoid the replace call when proper Expression support is in place
    #       and element/domain assignment is removed from compute_form_data.
    integrands = {
        num_points: replace(sorted_integrals[num_points].integrand(),
                            form_data.function_replace_map)
        for num_points in sorted(sorted_integrals)
    }

    # Build the more uflacs-specific intermediate representation
    uflacs_ir = build_uflacs_ir(itg_data.domain.ufl_cell(), itg_data.integral_type,
                                ir["entitytype"], integrands, ir["tensor_shape"],
                                quadrature_rules, parameters)

    # Add coefficient numbering to UFLACS IR
    uflacs_ir["coefficient_numbering"] = coefficient_numbering

    # Get coefficient offset
    coefficients = []
    factorization = uflacs_ir["piecewise_ir"]["F"]
    for i, attr in factorization.nodes.items():
        if attr['status'] == 'piecewise':
            v = attr['expression']
            mt = attr.get('mt', False)
            if mt and not v._ufl_is_literal_:
                if isinstance(mt.terminal, ufl.Coefficient):
                    tr = attr.get('tr', False)
                    begin, end = tr.dofrange
                    assert end - begin > 0
                    n = coefficient_numbering[mt.terminal]
                    # print("C0:", mt.terminal, n, end - begin)
                    coefficients.append((n, end - begin, mt.terminal))

    for num_points in uflacs_ir["all_num_points"]:
        factorization = uflacs_ir["varying_irs"][num_points]["F"]
        for i, attr in factorization.nodes.items():
            if attr['status'] == 'varying':
                v = attr['expression']
                mt = attr.get('mt', False)
                if mt and not v._ufl_is_literal_:
                    if isinstance(mt.terminal, ufl.Coefficient):
                        tr = attr.get('tr', False)
                        begin, end = tr.dofrange
                        assert end - begin > 0
                        n = coefficient_numbering[mt.terminal]
                        # print("C1:", mt.terminal, n, end - begin, mt.terminal.count())
                        coefficients.append((n, end - begin, mt.terminal))

    offsets = {}
    _offset = 0
    for num, size, c in sorted(coefficients):
        print(num, size, c)
        if c not in offsets:
            offsets[c] = _offset
            _offset += size
    print('Offsets = ', offsets.values())
    print('Numbering = ', coefficient_numbering.values())

    # Copy offsets into IR
    uflacs_ir["coefficient_offsets"] = offsets

    ir.update(uflacs_ir)

    return ir
