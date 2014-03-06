
import os
from uflacs.utils.log import info, warning

def add_compile_toy_options(opts):
    "Args: list of .ufl file(s)."
    pass

def run_compile_toy(options, args):
    "Compile forms and expressions from .ufl file(s) into C++ code with undefined FEM backend."
    from ufl.algorithms import load_ufl_file
    from uflacs.backends.toy.toy_compiler import compile_element, compile_expression, compile_form
    from uflacs.codeutils.format_code import format_code

    for input_filename in args:
        prefix, ext = os.path.splitext(os.path.basename(input_filename))
        if ext != '.ufl':
            warning("Expecting ufl file, got %s." % ext)
        output_filename = prefix + '.h'

        info("Loading file '%s'..." % (input_filename,))
        data = load_ufl_file(input_filename)

        info("Compiling '%s'..." % prefix)
        code = format_code([
            ['// Elements'],
            [compile_element(element, prefix) for element in data.elements],
            ['// Expressions'],
            [compile_expression(expr, prefix) for expr in data.expressions],
            ['// Forms'],
            [compile_form(form, prefix) for form in data.forms],
            ])

        info("Writing code to '%s'..." % output_filename)
        with open(output_filename, "w") as f:
            f.write(code)

    return 0

