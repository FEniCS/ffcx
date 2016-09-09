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

#from uflacs.params import default_parameters
from uflacs.analysis.modified_terminals import analyse_modified_terminal
from uflacs.representation.compute_expr_ir import compute_expr_ir
from uflacs.elementtables.terminaltables import build_element_tables, optimize_element_tables


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
    # uflacs_ir["name"] = form_data.name
    # uflacs_ir["coefficient_names"] = form_data.coefficient_names
    # uflacs_ir["argument_names"] = form_data.argument_names
    # uflacs_ir["cell"] = form_data.integration_domains[0].ufl_cell()
    # uflacs_ir["function_replace_map"] = form_data.function_replace_map

    # Build coefficient numbering for UFC interface here, to avoid renumbering in UFL and application of replace mapping
    uflacs_ir["coefficient_numbering"] = {}
    #uflacs_ir["coefficient_element"] = {}
    #uflacs_ir["coefficient_domain"] = {}
    for i, f in enumerate(sorted_by_count(form_data.function_replace_map.keys())):
        g = form_data.function_replace_map[f]
        assert i == g.count()
        uflacs_ir["coefficient_numbering"][g] = i # USING THIS ONE BECAUSE WE'RE CALLING REPLACE BELOW
        #uflacs_ir["coefficient_numbering"][f] = i # If we make elements and domains well formed we can avoid replace below and use this line instead
        #uflacs_ir["coefficient_element"][f] = g.ufl_element()
        #uflacs_ir["coefficient_domain"][f] = g.ufl_domain()

    # Build ir for each num_points/integrand
    uflacs_ir["expr_ir"] = {}
    for num_points in sorted(integrals_dict.keys()):
        integral = integrals_dict[num_points]

        # Get integrand
        expr = integral.integrand()

        # Replace coefficients so they all have proper element and domain for what's to come
        # TODO: We can avoid this step when Expression is in place and
        #       element/domain assignment removed from compute_form_data.
        # TODO: Doesn't replace domain coefficient!!!
        #       Merge replace functionality into change_to_reference_grad to fix?
        #       When coordinate field coefficient is removed I guess this issue will disappear?
        expr = replace(expr, form_data.function_replace_map) # FIXME: Still need to apply this mapping.

        # Build the core uflacs ir of expressions
        expr_ir = compute_expr_ir(expr)
        uflacs_ir["expr_ir"][num_points] = expr_ir

    for num_points in sorted(integrals_dict.keys()):
        expr_ir = uflacs_ir["expr_ir"][num_points]

        # Build set of modified terminal ufl expressions
        V = expr_ir["V"]
        modified_terminals = [analyse_modified_terminal(V[i])
                              for i in expr_ir["modified_terminal_indices"]]

        # Analyse modified terminals and store data about them
        terminal_data = modified_terminals + expr_ir["modified_arguments"]

        # Build tables needed by all modified terminals
        # (currently build here means extract from ffc psi_tables)
        #print '\n'.join([str(mt.expr) for mt in terminal_data])
        tables, terminal_table_names = build_element_tables(psi_tables, num_points,
                                                            entitytype, terminal_data,
                                                            epsilon)

        # Optimize tables and get table name and dofrange for each modified terminal
        unique_tables, terminal_table_ranges = optimize_element_tables(tables, terminal_table_names, epsilon)
        expr_ir["unique_tables"] = unique_tables

        # Modify ranges for restricted form arguments (not geometry!)
        # FIXME: Should not coordinate dofs get the same offset?
        from ufl.classes import FormArgument
        for i, mt in enumerate(terminal_data):
            # TODO: Get the definition that - means added offset from somewhere
            if mt.restriction == "-" and isinstance(mt.terminal, FormArgument):
                # offset = number of dofs before table optimization
                offset = int(tables[terminal_table_names[i]].shape[-1])
                (unique_name, b, e) = terminal_table_ranges[i]
                terminal_table_ranges[i] = (unique_name, b + offset, e + offset)

        # Split into arguments and other terminals before storing in expr_ir
        # TODO: Some tables are associated with num_points, some are not
        #       (i.e. piecewise constant, averaged and x0)
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
