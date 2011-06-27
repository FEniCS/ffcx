
import os
from uflacs.utils.log import info, warning

from uflacs.codeutils.format_code import format_code, Namespace
from uflacs.codeutils.dolfin_compiler import compile_dolfin_expression

def add_compile_dolfin_options(opts):
    "Args: list of .ufl file(s)."
    pass

def write_file(input_filename, output_filename, code):
    if not output_filename:
        output_filename = input_filename.replace('.ufl', '.h')
    f = open(output_filename, 'w')
    f.write(code)
    f.close()

def run_compile_dolfin(options, args):
    "Compile expressions from .ufl file(s) into dolfin C++ Expressions."
    from ufl.algorithms import load_ufl_file

    # TODO: Arguments?
    output_filename = None
    namespace = 'uflacs'
    includes = ['#include <iostream>', '#include <cmath>', '#include <dolfin.h>']

    # Generate code for each ufl file
    for filename in args:
        info("Loading file '%s'..." % (filename,))
        data = load_ufl_file(filename)

        prefix, ext = os.path.splitext(os.path.basename(filename))
        if ext != '.ufl':
            warning("Expecting ufl file, got %s." % ext)

        # Generate code for each expression in this file
        file_code = []
        for k, expr in enumerate(data.expressions):
            name = data.object_names.get(id(expr), 'w%d' % k)
            expr_code = compile_dolfin_expression(expr, name, data.object_names)
            file_code.append(expr_code)

        # Wrap code from each file in its own namespace
        code = format_code([includes, Namespace(prefix, file_code)])

        # Write it to file!
        write_file(filename, output_filename, code)

    return 0
