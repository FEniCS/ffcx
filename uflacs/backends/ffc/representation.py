# -*- coding: utf-8 -*-
# Copyright (C) 2011-2015 Martin Sandve Alnæs
#
# This file is part of UFLACS.
#
# UFLACS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# UFLACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with UFLACS. If not, see <http://www.gnu.org/licenses/>.

"""The FFC specific backend to the UFLACS form compiler algorithms."""

import numpy

from ufl.algorithms import replace
from ufl.utils.sorting import sorted_by_count
from ufl.classes import FormArgument, CellCoordinate

#from uflacs.params import default_parameters
from uflacs.analysis.modified_terminals import analyse_modified_terminal
from uflacs.representation.compute_expr_ir import compute_expr_ir
from uflacs.elementtables.terminaltables import build_element_tables, optimize_element_tables
from uflacs.backends.ffc.common import ufc_restriction_offset


def compute_uflacs_integral_ir(psi_tables, entitytype,
                               integrals_dict, form_data,
                               parameters):
    # TODO: Hack before we get default parameters properly into ffc
    #p = default_parameters()
    #p.update(parameters)
    #parameters = p

    # FIXME: Should be epsilon from ffc parameters
    from uflacs.language.format_value import get_float_threshold
    epsilon = get_float_threshold()

    uflacs_ir = {}

    # Some form_data info that we may need but currently don't use
    #form_data.name
    #form_data.coefficient_names
    #form_data.argument_names
    #form_data.integration_domains[0].ufl_cell()
    #form_data.function_replace_map

    # Get integrands (usually just one)
    all_num_points = sorted(integrals_dict.keys())
    integrands = {
        num_points: integrals_dict[num_points].integrand()
        for num_points in all_num_points
        }

    # Build coefficient numbering for UFC interface here, to avoid
    # renumbering in UFL and application of replace mapping
    sorted_coefficients = sorted_by_count(form_data.function_replace_map.keys())


    uflacs_ir["coefficient_numbering"] = {}
    if 0:
        pass
        # If we make elements and domains well formed we can
        # avoid replace below and use this code instead
        #uflacs_ir["coefficient_element"] = {}
        #uflacs_ir["coefficient_domain"] = {}
        #for i, f in enumerate(sorted_coefficients):
        #    g = form_data.function_replace_map[f]
        #    uflacs_ir["coefficient_numbering"][f] = i
        #    uflacs_ir["coefficient_element"][f] = g.ufl_element()
        #    uflacs_ir["coefficient_domain"][f] = g.ufl_domain()
    else:
        # Using this version because we're calling replace below
        for i, f in enumerate(sorted_coefficients):
            g = form_data.function_replace_map[f]
            uflacs_ir["coefficient_numbering"][g] = i
            assert i == g.count()
        # Replace coefficients so they all have proper element and domain for what's to come
        # TODO: We can avoid this step when proper Expression support is in place
        #       and element/domain assignment is removed from compute_form_data.
        integrands = {
            num_points: replace(integrands[num_points], form_data.function_replace_map)
            for num_points in all_num_points
            }


    # Build the core uflacs expression ir for each num_points/integrand
    # TODO: Better to compute joint IR for all integrands
    #       and deal with num_points later? If we want to
    #       adjoint quadrature rules for subterms automatically
    #       anyway, num_points should be advisory.
    #       For now, expecting multiple num_points to be rare.
    uflacs_ir["expr_irs"] = {
        num_points: compute_expr_ir(integrands[num_points])
        for num_points in all_num_points
        }


    for num_points in all_num_points:
        expr_ir = uflacs_ir["expr_irs"][num_points]

        # Build set of modified terminal ufl expressions
        V = expr_ir["V"]
        modified_terminals = [analyse_modified_terminal(V[i])
                              for i in expr_ir["modified_terminal_indices"]]
        terminal_data = modified_terminals + expr_ir["modified_arguments"]

        # Figure out if we need to access CellCoordinate to
        # avoid generating quadrature point table otherwise
        expr_ir["need_points"] = any(isinstance(mt.terminal, CellCoordinate)
                                     for mt in modified_terminals)

        # Build tables needed by all modified terminals
        # (currently build here means extract from ffc psi_tables)
        tables, terminal_table_names = \
            build_element_tables(psi_tables, num_points, entitytype, terminal_data, epsilon)

        # Optimize tables and get table name and dofrange for each modified terminal
        unique_tables, terminal_table_ranges = \
            optimize_element_tables(tables, terminal_table_names, epsilon)

        # Modify ranges for restricted form arguments
        # (geometry gets padded variable names instead)
        for i, mt in enumerate(terminal_data):
            if mt.restriction and isinstance(mt.terminal, FormArgument):
                # offset = 0 or number of dofs before table optimization
                num_original_dofs = int(tables[terminal_table_names[i]].shape[-1])
                offset = ufc_restriction_offset(mt.restriction, num_original_dofs)
                (unique_name, b, e) = terminal_table_ranges[i]
                terminal_table_ranges[i] = (unique_name, b + offset, e + offset)

        # Store the tables
        expr_ir["unique_tables"] = unique_tables

        # Split into arguments and other terminals before storing in expr_ir
        # TODO: Some tables are associated with num_points, some are not
        #       (i.e. piecewise constant, averaged and x0).
        #       It will be easier to deal with that if we can join
        #       the expr_ir for all num_points as mentioned above.
        n = len(expr_ir["modified_terminal_indices"])
        m = len(expr_ir["modified_arguments"])
        assert len(terminal_data) == n + m
        assert len(terminal_table_ranges) == n + m
        assert len(terminal_table_names) == n + m
        expr_ir["modified_terminal_table_ranges"] = terminal_table_ranges[:n]
        expr_ir["modified_argument_table_ranges"] = terminal_table_ranges[n:]

        # Store table data in V indexing, this is used in integralgenerator
        expr_ir["table_ranges"] = numpy.empty(len(V), dtype=object)
        expr_ir["table_ranges"][expr_ir["modified_terminal_indices"]] = \
            expr_ir["modified_terminal_table_ranges"]

    return uflacs_ir
