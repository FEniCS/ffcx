
from log import info
from utils import format_list, format_dict

def add_compile_options(opts):
    "Args: list of .ufl file(s)."
    pass

def compile_expression(expr):
    return "// TODO: format as C code: " + str(expr)

def compile_form(form_data):
    return "// TODO: format form as C code"

def run_compile(options, args):
    "Compile forms and expressions from .ufl file(s) into C expressions."
    from ufl.algorithms import load_ufl_file
    filenames = args
    code = []
    for fn in filenames:
        info("Loading file '%s'..." % (fn,))
        data = load_ufl_file(fn)

        code.append([\
            ['// Expressions'],
            [compile_expression(expr) for expr in data.expressions],
            ['// Forms'],
            [compile_form(form) for form in data.forms],
            ])

    print '\n\n'.join([format_list(c) for c in code])

    return 0
