
from uflacs.utils.log import uflacs_assert
from uflacs.analysis.dependency_handler import DependencyHandler
from uflacs.codeutils.format_code_structure import format_code_structure, Indented, ArrayDecl
from uflacs.generation.generate import generate_code_from_ssa, generate_expression_body
from uflacs.backends.ffc.ffc_language_formatter import FFCLanguageFormatter
from uflacs.backends.ffc.ffc_statement_formatter import FFCStatementFormatter

def generate_tabulate_tensor_code(ir, parameters):

    # Fetch uflacs specific ir part
    uflacs_ir = ir["uflacs"]

    # Create an object to track dependencies across other components
    dependency_handler = DependencyHandler(uflacs_ir["terminals"],
                                           uflacs_ir["function_replace_map"],
                                           uflacs_ir["argument_names"],
                                           uflacs_ir["coefficient_names"])

    # Create backend specific plugin objects
    language_formatter = FFCLanguageFormatter(dependency_handler, ir)
    statement_formatter = FFCStatementFormatter(dependency_handler, ir)

    # Generate code partitions from ir
    all_num_points = sorted(uflacs_ir["expression_partitions"])
    for num_points in all_num_points:
        partitions_ir = uflacs_ir["expression_partitions"][num_points]

        # TODO: Possibly need to pass the full uflacs_ir["expression_partitions"] dict to generate_code_...
        uflacs_assert(len(all_num_points) == 1,
                      "Assuming a single quadrature rule per integral domain from this point on in uflacs.")

        partition_codes, final_variable_names = generate_code_from_ssa(partitions_ir, language_formatter)

    # Generate full code from snippets
    expression_body = generate_expression_body(statement_formatter,
                                               partition_codes,
                                               final_variable_names,
                                               uflacs_ir["num_registers"])

    # Format uflacs specific code structures into a single
    # string and place in dict before returning to ffc
    body = format_code_structure(Indented([
        language_formatter.get_using_statements(),
        uflacs_ir["tables_code"], # FIXME: Refactoring of build_tables needed to generate table code here
        expression_body,
        ]))
    code = {
        "tabulate_tensor": body,
        "additional_includes_set": language_formatter.get_includes(),
        }
    return code


def new_generate_tabulate_tensor_code(ir, parameters):

    # Fetch uflacs specific ir part
    uflacs_ir = ir["uflacs"]

    # Create backend specific plugin objects
    language_formatter = FFCLanguageFormatter(dependency_handler, ir)
    statement_formatter = FFCStatementFormatter(dependency_handler, ir)

    # Generate code partitions from ir
    all_num_points = sorted(uflacs_ir["expression_partitions"])
    for num_points in all_num_points:
        partitions_ir = uflacs_ir["expression_partitions"][num_points]

        # TODO: Possibly need to pass the full uflacs_ir["expression_partitions"] dict to generate_code_...
        uflacs_assert(len(all_num_points) == 1,
                      "Assuming a single quadrature rule per integral domain from this point on in uflacs.")

        partition_codes, final_variable_names = generate_code_from_ssa(partitions_ir, language_formatter)

    # Generate full code from snippets
    expression_body = generate_expression_body(statement_formatter,
                                               partition_codes,
                                               final_variable_names,
                                               uflacs_ir["num_registers"])

    # Format uflacs specific code structures into a single
    # string and place in dict before returning to ffc
    body = format_code_structure(Indented([
        language_formatter.get_using_statements(),
        uflacs_ir["tables_code"], # FIXME: Refactoring of build_tables needed to generate table code here
        expression_body,
        ]))
    code = {
        "tabulate_tensor": body,
        "additional_includes_set": language_formatter.get_includes(),
        }
    return code
