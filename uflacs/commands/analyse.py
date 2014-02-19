
from uflacs.utils.log import info
from uflacs.utils.str_utils import format_list

def add_analyse_options(opts):
    "Args: list of .ufl file(s)."
    pass

def run_analyse(options, args):
    "Analyse various properties of forms found in .ufl file(s)."
    import ufl
    from ufl.algorithms import load_ufl_file
    ufl.algorithms.preprocess.enable_profiling = True
    filenames = args
    for fn in filenames:
        info("Loading file '%s'..." % (fn,))
        data = load_ufl_file(fn)
        for form in data.forms:
            form_data = form.compute_form_data()

            # Stuff we can analyse: TODO: Go through list and implement
            # cell
            # elements
            # geometric quantities
            # coefficients
            # expressions
            #   shape
            #     argument deps
            #     geometry deps
            #     coefficient deps
            #     polynomial degree
            # forms
            #   integrals
            #     measure
            #     arguments
            #     geometry deps
            #     coefficient deps
            #     polynomial degree

            info("Families: \n" + format_list([e.family() for e in form_data.elements]))
            info("Elements: \n" + format_list(form_data.elements))
    return 0
