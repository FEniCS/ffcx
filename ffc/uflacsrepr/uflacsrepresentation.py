# Copyright (C) 2013-2014 Martin Alnaes
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

from ffc.log import info, error, begin, end, debug_ir, ffc_assert, warning

from ffc.fiatinterface import create_element
from ffc.representationutils import initialize_integral_ir
from ffc.quadrature.parameters import parse_optimise_parameters
from ffc.quadrature.tabulate_basis import tabulate_basis
from ffc.quadrature.quadraturerepresentation import sort_integrals

from ufl.algorithms import replace
from ufl.algorithms.change_to_reference import (change_to_reference_grad,
                                                compute_integrand_scaling_factor,
                                                change_to_reference_geometry)
from ufl.algorithms import propagate_restrictions

from uflacs.utils.log import uflacs_assert
from uflacs.params import default_parameters
from uflacs.analysis.datastructures import object_array
from uflacs.analysis.modified_terminals import analyse_modified_terminal2
from uflacs.representation.compute_expr_ir import compute_expr_ir
from uflacs.elementtables.terminaltables import build_element_tables, optimize_element_tables

def compute_integral_ir(itg_data,
                        form_data,
                        form_id,
                        parameters):
    "Compute intermediate represention of integral."

    info("Computing uflacs representation")

    # Initialise representation
    ir = initialize_integral_ir("uflacs", itg_data, form_data, form_id)

    # Sort integrals into a dict with quadrature degree and rule as key
    sorted_integrals = sort_integrals(itg_data.integrals,
                                      itg_data.metadata["quadrature_degree"],
                                      itg_data.metadata["quadrature_rule"])

    # TODO: Might want to create the uflacs ir first and then create the tables we need afterwards!
    # Tabulate quadrature points and basis function values in these points
    integrals_dict, psi_tables, quadrature_rules = \
        tabulate_basis(sorted_integrals, form_data, itg_data)

    # Delegate to flacs to build its intermediate representation and add to ir
    uflacs_ir = compute_tabulate_tensor_ir(psi_tables, ir["entitytype"], integrals_dict, form_data, parameters)

    # Store uflacs generated part separately
    ir["uflacs"] = uflacs_ir

    # Create and save the optisation parameters # TODO: Define uflacs specific optimization parameters instead
    #ir["optimise_parameters"] = parse_optimise_parameters(parameters)

    # Save tables for quadrature weights and points
    ir["quadrature_rules"] = quadrature_rules

    # Create dimensions of primary indices, needed to reset the argument 'A'
    # given to tabulate_tensor() by the assembler.
    ir["prim_idims"] = [create_element(ufl_element).space_dimension()
                        for ufl_element in form_data.argument_elements]

    # Added for uflacs, not sure if this is the best way to get this:
    ir["coeff_idims"] = [create_element(ufl_element).space_dimension()
                         for ufl_element in form_data.coefficient_elements]

    return ir

def compute_tabulate_tensor_ir(psi_tables, entitytype,
                               integrals_dict, form_data,
                               parameters):
    # TODO: Hack before we get default parameters properly into ffc
    p = default_parameters()
    p.update(parameters)
    parameters = p

    uflacs_ir = {}
    #uflacs_ir["name"] = form_data.name
    #uflacs_ir["coefficient_names"] = form_data.coefficient_names
    #uflacs_ir["argument_names"] = form_data.argument_names
    #uflacs_ir["cell"] = form_data.integration_domains[0].cell()
    #uflacs_ir["function_replace_map"] = form_data.function_replace_map

    # Build ir for each num_points/integrand
    uflacs_ir["expr_ir"] = {}
    for num_points in sorted(integrals_dict.keys()):
        integral = integrals_dict[num_points]


        # FIXME: Move this symbolic processing to compute_form_data, give compute_form_data an option to to this for now.


        # Get integrand expr and apply some symbolic preprocessing
        expr = integral.integrand()

        # Replace coefficients so they all have proper element and domain for what's to come
        expr = replace(expr, form_data.function_replace_map)

        # Change from physical gradients to reference gradients
        expr = change_to_reference_grad(expr) # TODO: Make this optional depending on backend

        # Compute and apply integration scaling factor
        scale = compute_integrand_scaling_factor(integral.domain(), integral.integral_type())
        expr = expr * scale

        # Change geometric representation to lower level quantities
        if integral.integral_type() == "quadrature":
            physical_coordinates_known = True
        else:
            physical_coordinates_known = False
        expr = change_to_reference_geometry(expr, physical_coordinates_known)

        # Restrictions may not be at terminals any more, propagate them again TODO: Skip this in preprocess?
        if integral.integral_type() == "interior_facet":
            expr = propagate_restrictions(expr)



        # Build the core uflacs ir of expressions
        expr_ir = compute_expr_ir(expr, parameters)
        uflacs_ir["expr_ir"][num_points] = expr_ir

    # NB! Using the last num_points from integrals_dict below, but not handling it properly yet
    uflacs_assert(len(integrals_dict) == 1, "Not supporting multiple integration rules yet.")
    num_points, = list(integrals_dict.keys())
    expr_ir = uflacs_ir["expr_ir"][num_points]

    # Build set of modified terminal ufl expressions
    V = expr_ir["V"]
    modified_terminals = [V[i] for i in expr_ir["modified_terminal_indices"]]
    modified_terminals += expr_ir["modified_arguments"]

    # Analyse modified terminals and store data about them
    terminal_data = [analyse_modified_terminal2(o) for o in modified_terminals]

    # Build tables needed by all modified terminals (currently build here means extract from ffc psi_tables)
    tables, terminal_table_names = build_element_tables(psi_tables, num_points,
                                                        entitytype, terminal_data)

    # Optimize tables and get table name and dofrange for each modified terminal
    unique_tables, terminal_table_ranges = optimize_element_tables(tables, terminal_table_names)
    expr_ir["unique_tables"] = unique_tables

    # Split into arguments and other terminals before storing in expr_ir
    # TODO: Some tables are associated with num_points, some are not (i.e. piecewise constant, averaged and x0)
    n = len(expr_ir["modified_terminal_indices"])
    m = len(expr_ir["modified_arguments"])
    assert len(terminal_table_ranges) == n+m
    assert len(terminal_table_names) == n+m
    expr_ir["modified_terminal_table_ranges"] = terminal_table_ranges[:n]
    expr_ir["modified_argument_table_ranges"] = terminal_table_ranges[n:]

    # Store table data in V indexing, this is used in integralgenerator
    expr_ir["table_ranges"] = object_array(len(V))
    expr_ir["table_ranges"][expr_ir["modified_terminal_indices"]] = expr_ir["modified_terminal_table_ranges"]
    #for i in xrange(n):
    #    expr_ir["table_ranges"][expr_ir["modified_terminal_indices"][i]] = terminal_table_ranges[i]

    assert len(expr_ir["modified_argument_table_ranges"]) == len(expr_ir["modified_arguments"])

    return uflacs_ir
