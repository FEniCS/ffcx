
from uflacs.utils.log import info
from uflacs.utils.str_utils import format_list, format_dict

def add_analyse_options(opts):
    "Args: list of .ufl file(s)."
    pass

def run_analyse(options, args):
    "Analyse various properties of forms found in .ufl file(s)."
    from ufl.algorithms import load_ufl_file
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


def add_testelementmapping_options(opts):
    pass

def run_testelementmapping(options, args):
    from ufl import triangle, quadrilateral
    from ufl.algorithms import load_ufl_file, extract_elements
    filenames = args
    for fn in filenames:
        print "Loading file '%s'..." % (fn,)
        data = load_ufl_file(fn)
        for form in data.forms:
            # Build a hack of an element mapping for testing
            element_mapping = {}
            for e in extract_elements(form):
                e2 = e.reconstruct(family="Discontinuous Lagrange", cell=quadrilateral)
                element_mapping[e] = e2
            print element_mapping

            form.compute_form_data(element_mapping=element_mapping)
            form_data = form.form_data()
            print form_data.elements
            print [e.family() for e in form_data.elements]
    return 0
