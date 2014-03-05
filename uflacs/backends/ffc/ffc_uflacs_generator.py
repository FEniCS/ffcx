
from uflacs.generation.integralgenerator import IntegralGenerator
from uflacs.backends.ffc.scratch import FFCLanguageFormatter

def generate_tabulate_tensor_code(ir, parameters):

    # Create FFC backend
    language_formatter = FFCLanguageFormatter(ir)

    # Create code generator for integral body
    ig = IntegralGenerator(ir, language_formatter)

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
