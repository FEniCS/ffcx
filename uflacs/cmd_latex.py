
def add_latex_options(opts):
    pass

def run_latex(options, args):
    print "Running latex"
    from uflacs.codeutils.latex_format_test import test_latex_formatting
    return test_latex_formatting()
