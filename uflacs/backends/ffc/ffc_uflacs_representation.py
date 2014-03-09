
from ufl.algorithms import replace
from ufl.algorithms.change_to_reference import (change_to_reference_grad,
                                                compute_integrand_scaling_factor,
                                                change_to_reference_geometry)

from uflacs.utils.log import uflacs_assert
from uflacs.params import default_parameters
from uflacs.analysis.datastructures import object_array
from uflacs.analysis.modified_terminals import analyse_modified_terminal2
from uflacs.representation.compute_expr_ir import compute_expr_ir
from uflacs.elementtables.terminaltables import build_element_tables, optimize_element_tables

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

        # Get integrand expr and apply some symbolic preprocessing
        expr = integral.integrand()

        # Replace coefficients so they all have proper element and domain for what's to come
        expr = replace(expr, form_data.function_replace_map)

        # Change from physical gradients to reference gradients
        expr = change_to_reference_grad(expr) # TODO: Make this optional depending on backend

        # Compute and apply integration scaling factor
        scale = compute_integrand_scaling_factor(integral.domain(), integral.domain_type())
        expr = expr * scale

        # Change geometric representation to lower level quantities
        if integral.domain_type() == "quadrature":
            physical_coordinates_known = True
        else:
            physical_coordinates_known = False
        expr = change_to_reference_geometry(expr, physical_coordinates_known)

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


def optimize_tabulate_tensor_ir(ir, parameters):
    # Hack before we get default parameters properly into ffc
    p = default_parameters()
    p.update(parameters)
    parameters = p

    # TODO: Implement some optimization here!
    oir = ir
    return oir
