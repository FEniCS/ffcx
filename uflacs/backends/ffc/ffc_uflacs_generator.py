
from uflacs.utils.log import uflacs_assert
from uflacs.analysis.dependency_handler import DependencyHandler
from uflacs.codeutils.format_code_structure import format_code_structure, Indented, ArrayDecl
from uflacs.generation.generate import generate_code_from_ssa, generate_expression_body
from uflacs.backends.ffc.ffc_language_formatter import FFCLanguageFormatter
from uflacs.backends.ffc.ffc_statement_formatter import FFCStatementFormatter

def generate_tabulate_tensor_code(ir, parameters):

    # FIXME: Create FFC backend
    #ffc_backend = FFCBackend(ir, parameters)

    # Fetch uflacs specific ir part
    ig = IntegralGenerator(ir) #, ffc_backend)

    # Generate code for the tabulate_tensor body
    body = ig.generate()

    # Fetch includes
    includes = [] # sorted(set(ig.get_includes() + backend.get_includes())) # FIXME

    # Format uflacs specific code structures into a single
    # string and place in dict before returning to ffc
    code = {
        "tabulate_tensor": body,
        "additional_includes_set": includes,
        }
    return code
