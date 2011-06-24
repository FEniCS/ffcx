
from uflacs.utils.log import info
from uflacs.utils.str_utils import format_list, format_dict

def add_latex_options(opts):
    "Args: list of .ufl file(s)."
    pass

def compile_element(element):
    return "TODO: format as LaTeX code: " + str(element)

def compile_expression(expr):
    return "TODO: format as LaTeX code: " + str(expr)

from uflacs.codeutils.latex_compiler import compile_form

def run_latex(options, args):
    "Compile forms and expressions from .ufl file(s) into LaTeX expressions."
    from uflacs.codeutils.format_code import format_code
    from ufl.algorithms import load_ufl_file
    filenames = args
    code = []
    for fn in filenames:
        info("Loading file '%s'..." % (fn,))
        data = load_ufl_file(fn)

        code.append([\
            [r'\section{Elements}'],
            [compile_element(element) for element in data.elements],
            [r'\section{Expressions}'],
            [compile_expression(expr) for expr in data.expressions],
            [r'\section{Forms}'],
            [compile_form(form) for form in data.forms],
            ])

    print format_code(code)

    return 0
