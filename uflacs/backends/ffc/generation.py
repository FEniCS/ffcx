
from uflacs.codeutils.cpp_expr_formatting_rules import CppExprFormatter
from uflacs.generation.integralgenerator import IntegralGenerator
from uflacs.backends.ffc.access import FFCAccessBackend
from uflacs.backends.ffc.definitions import FFCDefinitionsBackend


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
    includes = set()
    includes.update(ig.get_includes())
    includes.update(backend_definitions.get_includes())

    # Format uflacs specific code structures into a single
    # string and place in dict before returning to ffc
    code = {
        "tabulate_tensor": body,
        "additional_includes_set": includes,
    }
    return code
