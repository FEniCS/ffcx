
import os
from uflacs.utils.log import info, warning

def add_compile_dolfin_options(opts):
    "Args: list of .ufl file(s)."
    pass # TODO: Output filename argument

def write_file(input_filename, output_filename, code):
    if not output_filename:
        output_filename = input_filename.replace('.ufl', '.h')
    f = open(output_filename, 'w')
    f.write(code)
    f.close()

def run_compile_dolfin(options, args):
    "Compile expressions from .ufl file(s) into dolfin C++ Expressions."
    from ufl.algorithms import load_ufl_file
    from uflacs.codeutils.dolfin_compiler import compile_dolfin_expressions_header

    # TODO: Arguments?
    output_filename = None

    # Generate code for each ufl file
    for filename in args:
        info("Loading file '%s'..." % (filename,))
        prefix, ext = os.path.splitext(os.path.basename(filename))
        if ext != '.ufl':
            warning("Expecting ufl file, got %s." % ext)

        data = load_ufl_file(filename)

        code = compile_dolfin_expressions_header(data, prefix)

        write_file(filename, output_filename, code)

    return 0
