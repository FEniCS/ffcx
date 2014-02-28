
from ufl.common import product
from ufl.classes import Argument, FormArgument
from ufl.algorithms import extract_unique_elements
from uflacs.utils.log import uflacs_assert
from uflacs.params import default_parameters
from uflacs.analysis.modified_terminals import analyse_modified_terminal #, analyse_modified_terminal2
from uflacs.analysis.table_utils import (generate_psi_table_name,
                                         get_ffc_table_values,
                                         strip_table_zeros,
                                         build_unique_tables,
                                         derivative_listing_to_counts,
                                         flatten_component)

from uflacs.generation.compiler import compile_expression_partitions

from uflacs.codeutils.format_code_structure import ArrayDecl

def build_element_counter_map(integrals_dict, element_replace_map):
    element_counter_map = {}
    for num_points in sorted(integrals_dict.keys()):
        element_counter_map[num_points] = {}
        ecm = element_counter_map[num_points]

        # Find all elements in this integrand and map them
        integrand = integrals_dict[num_points].integrand()
        elements = [element_replace_map[e] for e in extract_unique_elements(integrand)]

        # Count the elements in a stable sorting
        for element in sorted(elements):
            if element not in ecm:
                ecm[element] = len(ecm)
    return element_counter_map

def build_tables(psi_tables, entitytype, element_counter_map, terminal_data):
    num_points, = element_counter_map.keys() # Assuming a single num_points value
    tables = {}
    preserve_tables = set()
    argument_tables = {}
    handled = set()
    for t, c, d, r in terminal_data:

        # Avoid duplicating tables because of restriction
        key = (t, c, d)
        if key in handled:
            continue
        handled.add(key)

        if isinstance(t, FormArgument):
            element = t.element()

            #import pdb
            #pdb.set_trace()

            element_counter = element_counter_map[num_points][element]

            # Change derivatives format for table lookup
            gdim = element.cell().geometric_dimension()
            derivatives = tuple(derivative_listing_to_counts(d, gdim))

            # Flatten component
            flat_component = flatten_component(c, t.shape(), element.symmetry())

            # Get name and values for this particular table
            name = generate_psi_table_name(element_counter, flat_component,
                                         derivatives, entitytype)
            table = get_ffc_table_values(psi_tables, entitytype, num_points,
                                         element, flat_component, derivatives)
            tables[name] = table
            if 0:
                print
                print element
                print name
                print table
                print

        if isinstance(t, Argument):
            # Avoid deleting the original table so we can loop over all dofs for arguments:
            preserve_tables.add(name)
            # Group argument tables by t,c for nonzero detection TODO: Not used yet
            if (t,c) not in argument_tables:
                argument_tables[(t,c)] = set()
            argument_tables[(t,c)].add(name)

    # FIXME: Build argument component dof ranges, here or somewhere else?
    #element_counter = element_counter_map[num_points][element]
    #gdim = element.cell().geometric_dimension()
    #derivatives = tuple(derivative_listing_to_counts(d, gdim))
    #flat_component = flatten_component(c, t.shape(), element.symmetry())
    #name = generate_psi_table_name(element_counter, flat_component,
    #                               derivatives, entitytype)
    #(uname, begin, end) = ir["table_ranges"][name]


    # Apply zero stripping to all tables
    stripped_tables = {}
    table_ranges = {}
    for name, table in tables.iteritems():
        begin, end, stable = strip_table_zeros(table)
        stripped_tables[name] = stable
        table_ranges[name] = (begin, end)

        # Preserve some tables under a modified name:
        if name in preserve_tables:
            pname = "p" + name # Hack! TODO: Make a cleaner solution!
            begin, end = 0, table.shape[-1]
            stripped_tables[pname] = table
            table_ranges[pname] = (begin, end)

    # Build unique table mapping
    unique_tables, table_mapping = build_unique_tables(stripped_tables)

    # Build mapping of constructed table names to unique names
    unique_table_names = {}
    mapping_to_name = {}
    for name in sorted(table_mapping.keys()):
        ui = table_mapping[name]
        if ui not in mapping_to_name:
            mapping_to_name[ui] = name
        uname = mapping_to_name[ui]
        unique_table_names[name] = uname

    # Build mapping of constructed table names to data: name -> unique name, table range begin, table range end
    table_data = {}
    for name in sorted(table_mapping.keys()):
        uname = unique_table_names[name]
        b, e = table_ranges[name]
        table_data[name] = (uname, b, e)

    # FIXME: Refactor to extract this table code generation!
    # Format unique tables into code
    tables_code = [ArrayDecl("static const double", mapping_to_name[ui],
                             table.shape, table)
                   for ui, table in enumerate(unique_tables)
                   if product(table.shape) > 0]
    return tables_code, table_data

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

    # Build num_points/element to counter mapping
    uflacs_ir["element_map"] = build_element_counter_map(integrals_dict, form_data.element_replace_map)

    # Build ir for each num_points/integrand
    uflacs_ir["expression_partitions"] = {}
    uflacs_assert(len(integrals_dict) == 1, "Not supporting multiple integration rules yet.")
    for num_points in sorted(integrals_dict.keys()):
        integrand = integrals_dict[num_points].integrand()

        # Build the main uflacs ir
        partitions_ir = compile_expression_partitions(integrand, form_data.function_replace_map, parameters)
        uflacs_ir["expression_partitions"][num_points] = partitions_ir

        modified_terminals = partitions_ir["terminals"] # set of modified terminal ufl expression

        # Analyse modified terminals and store data about them in a canonical ordering
        terminal_data = [analyse_modified_terminal(o, uflacs_ir["function_replace_map"])
                         for o in sorted_expr(modified_terminals)]


    # Analyse the psi_tables that are required by functions etc.
    tables_code, table_data = build_tables(ir["psi_tables"],
                                           ir["entitytype"],
                                           uflacs_ir["element_map"],
                                           terminal_data)
    uflacs_ir["table_ranges"] = table_data
    uflacs_ir["tables_code"] = tables_code

    return uflacs_ir


def optimize_tabulate_tensor_ir(ir, parameters):
    # Hack before we get default parameters properly into ffc
    p = default_parameters()
    p.update(parameters)
    parameters = p

    # TODO: Implement some optimization here!
    oir = ir
    return oir
