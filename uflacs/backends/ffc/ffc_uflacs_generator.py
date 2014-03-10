
from uflacs.codeutils.cpp_expr_formatting_rules import CppExprFormatter
from uflacs.backends.ffc.ffc_backend import FFCAccessBackend, FFCDefinitionsBackend
from uflacs.generation.integralgenerator import IntegralGenerator

def generate_tabulate_tensor_code(ir, parameters):

    # Create C++ backend
    language_formatter = CppExprFormatter()

    # Create FFC backend
    backend_access = FFCAccessBackend(ir, parameters)
    backend_definitions = FFCDefinitionsBackend(ir, parameters)

    # Create code generator for integral body
    ig = IntegralGenerator(ir, language_formatter, backend_access, backend_definitions)

    # Generate code for the tabulate_tensor body
    body = ig.generate()

    # Fetch includes
    includes = [] # sorted(set(ig.get_includes() + backend.get_includes())) # FIXME
    includes = set(("#include <cstring>",
                    "#include <cmath>"))

    # Format uflacs specific code structures into a single
    # string and place in dict before returning to ffc
    code = {
        "tabulate_tensor": body,
        "additional_includes_set": includes,
        }
    return code
