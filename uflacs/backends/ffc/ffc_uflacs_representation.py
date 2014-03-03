
from uflacs.utils.log import uflacs_assert
from uflacs.params import default_parameters
from uflacs.analysis.modified_terminals import analyse_modified_terminal2
from uflacs.generation.compiler import compute_expr_ir

from ufl.algorithms import replace, change_to_local_grad
def compute_tabulate_tensor_ir(integrals_dict,
                               form_data,
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
        # Get integrand expr and apply some symbolic preprocessing
        expr = integrals_dict[num_points].integrand()
        expr = replace(expr, form_data.function_replace_map)
        expr = change_to_local_grad(expr) # TODO: Make this optional depending on backend

        # Build the core uflacs ir of expressions
        expr_ir = compute_expr_ir(expr, parameters)
        uflacs_ir["expr_ir"][num_points] = expr_ir

    # NB! Using the last num_points from integrals_dict below, but not handling it properly yet
    uflacs_assert(len(integrals_dict) == 1, "Not supporting multiple integration rules yet.")
    num_points, = list(integrals_dict.keys())
    expr_ir = ir["expr_ir"][num_points]

    # Build set of modified terminal ufl expressions
    modified_terminals = [V[i] for i in expr_ir["modified_terminal_indices"]]
    modified_terminals += expr_ir["modified_arguments"]

    # Analyse modified terminals and store data about them
    terminal_data = [analyse_modified_terminal2(o) for o in modified_terminals]

    # Build tables needed by all modified terminals (currently build here means extract from ffc psi_tables)
    tables, terminal_table_names = build_element_tables(ir["psi_tables"], num_points,
                                                        ir["entitytype"], terminal_data)

    # Optimize tables and get table name and dofrange for each modified terminal
    unique_tables, terminal_table_ranges = optimize_element_tables(tables, terminal_table_names)

    # Split into arguments and other terminals before storing in expr_ir
    # TODO: Some tables are associated with num_points, some are not (i.e. piecewise constant, averaged and x0)
    n = len(modified_terminals)
    expr_ir["modified_terminal_table_ranges"] = terminal_table_ranges[:n]
    expr_ir["modified_terminal_table_ranges"] = terminal_table_ranges[n:]

    return uflacs_ir


def optimize_tabulate_tensor_ir(ir, parameters):
    # Hack before we get default parameters properly into ffc
    p = default_parameters()
    p.update(parameters)
    parameters = p

    # TODO: Implement some optimization here!
    oir = ir
    return oir
