
import os
from ffc.log import info, warning

def add_latex_options(opts):
    "Args: list of .ufl file(s)."
    pass

def run_latex(options, args):
    "Compile forms and expressions from .ufl file(s) into LaTeX expressions."
    from ufl.algorithms import load_ufl_file
    from uflacs.backends.latex.compiler import compile_latex_document

    for filename in args:
        info("Loading file '%s'..." % (filename,))

        prefix, ext = os.path.splitext(os.path.basename(filename))
        if ext != '.ufl':
            warning("Expecting ufl file, got %s." % ext)
        output_filename = prefix + '.tex'

        data = load_ufl_file(filename)
        code = compile_latex_document(data)
        with open(output_filename, 'w') as f:
            f.write(code)

    return 0
