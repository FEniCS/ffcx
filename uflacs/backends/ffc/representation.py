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

from ufl.algorithms import propagate_restrictions
from ufl.algorithms import replace
from ufl.algorithms.change_to_reference import (change_to_reference_grad,
                                                compute_integrand_scaling_factor,
                                                change_to_reference_geometry)

from ffc.log import ffc_assert

from uflacs.params import default_parameters
from uflacs.datastructures.arrays import object_array
from uflacs.analysis.modified_terminals import analyse_modified_terminal
from uflacs.representation.compute_expr_ir import compute_expr_ir
from uflacs.elementtables.terminaltables import build_element_tables, optimize_element_tables


def compute_uflacs_integral_ir(psi_tables, entitytype,
                               integrals_dict, form_data,
                               parameters):
    # TODO: Hack before we get default parameters properly into ffc
    p = default_parameters()
    p.update(parameters)
    parameters = p

    uflacs_ir = {}
    # uflacs_ir["name"] = form_data.name
    # uflacs_ir["coefficient_names"] = form_data.coefficient_names
    # uflacs_ir["argument_names"] = form_data.argument_names
    # uflacs_ir["cell"] = form_data.integration_domains[0].cell()
    # uflacs_ir["function_replace_map"] = form_data.function_replace_map

    # Build ir for each num_points/integrand
    uflacs_ir["expr_ir"] = {}
    for num_points in sorted(integrals_dict.keys()):
        integral = integrals_dict[num_points]

        # FIXME: Move this symbolic processing to compute_form_data, give compute_form_data an option to to this for now.
        # Get integrand expr and apply some symbolic preprocessing
        expr = integral.integrand()

        # Replace coefficients so they all have proper element and domain for what's to come
        expr = replace(expr, form_data.function_replace_map)  # FIXME: Doesn't replace domain coefficient!!! Merge replace functionality into change_to_reference_grad to fix?

        # Change from physical gradients to reference gradients
        expr = change_to_reference_grad(expr)  # TODO: Make this optional depending on backend

        # Compute and apply integration scaling factor
        scale = compute_integrand_scaling_factor(integral.domain(), integral.integral_type())
        expr = expr * scale

        # Change geometric representation to lower level quantities
        if integral.integral_type() in ("custom", "vertex"):
            physical_coordinates_known = True
        else:
            physical_coordinates_known = False
        expr = change_to_reference_geometry(expr, physical_coordinates_known, form_data.function_replace_map)

        # Propagate restrictions to become terminal modifiers because change_to_reference_geometry
        # may have messed it up after the first pass in compute_form_data...
        if integral.integral_type() == "interior_facet":
            expr = propagate_restrictions(expr)

        # Build the core uflacs ir of expressions
        expr_ir = compute_expr_ir(expr, parameters)
        uflacs_ir["expr_ir"][num_points] = expr_ir

    # NB! Using the last num_points from integrals_dict below, but not handling it properly yet
    #ffc_assert(len(integrals_dict) == 1, "Not supporting multiple integration rules yet.")
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
        tables, terminal_table_names = build_element_tables(psi_tables, num_points,
                                                            entitytype, terminal_data)

        # Optimize tables and get table name and dofrange for each modified terminal
        unique_tables, terminal_table_ranges = optimize_element_tables(tables, terminal_table_names)
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
        expr_ir["table_ranges"] = object_array(len(V))
        expr_ir["table_ranges"][expr_ir["modified_terminal_indices"]] = \
            expr_ir["modified_terminal_table_ranges"]

    return uflacs_ir
