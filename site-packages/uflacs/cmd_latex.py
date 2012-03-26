
import os
from uflacs.utils.log import info, warning

def add_latex_options(opts):
    "Args: list of .ufl file(s)."
    pass # TODO: Output filename argument

def write_file(input_filename, output_filename, code):
    if not output_filename:
        output_filename = input_filename.replace('.ufl', '.tex')
    f = open(output_filename, 'w')
    f.write(code)
    f.close()

def run_latex(options, args):
    "Compile forms and expressions from .ufl file(s) into LaTeX expressions."
    from ufl.algorithms import load_ufl_file
    #from uflacs.codeutils.latex_compiler import compile_latex_document

    # TODO: Arguments?
    output_filename = None

    # Generate code for each ufl file
    for filename in args:
        info("Loading file '%s'..." % (filename,))
        prefix, ext = os.path.splitext(os.path.basename(filename))
        if ext != '.ufl':
            warning("Expecting ufl file, got %s." % ext)

        data = load_ufl_file(filename)

        code = compile_latex_document(data)

        write_file(filename, output_filename, code)

    return 0

# TODO: Move all these to latex compiler file
def compile_latex_document(data):
    from uflacs.codeutils.format_code import format_code
    from uflacs.codeutils.latex_compiler import compile_form
    code = [\
        [r'\section{Elements}'],
        [compile_element(element) for element in data.elements],
        [r'\section{Expressions}'],
        [compile_expression(expr) for expr in data.expressions],
        [r'\section{Forms}'],
        [compile_form(form) for form in data.forms],
        ]
    return format_code(code)

def compile_element(element):
    return "TODO: format as LaTeX code: " + str(element)

def compile_expression(expr):
    return "TODO: format as LaTeX code: " + str(expr)
