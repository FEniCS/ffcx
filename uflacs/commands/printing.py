
from ffc.log import info

def load_and_print_objects(filenames, tostr, skip_elements=False,
                           skip_coefficients=False, skip_expressions=False,
                           skip_forms=False):
    from ufl.algorithms import load_ufl_file
    for fn in filenames:
        info("Loading file '%s'..." % (fn,))
        data = load_ufl_file(fn)
        if not skip_elements:
            for o in data.elements:
                name = data.object_names.get(id(o), "<undefined element>")
                info("Element '%s':\n%s\n" % (name, tostr(o)))
        if not skip_coefficients:
            for o in data.coefficients:
                name = data.object_names.get(id(o), "<undefined coefficient>")
                info("Coefficient '%s':\n%s\n" % (name, tostr(o)))
        if not skip_expressions:
            for o in data.expressions:
                name = data.object_names.get(id(o), "<undefined expression>")
                info("Expr '%s':\n%s\n" % (name, tostr(o)))
        if not skip_forms:
            for o in data.forms:
                name = data.object_names.get(id(o), "<undefined form>")
                info("Form '%s':\n%s\n" % (name, tostr(o)))
    return 0


def add_str_options(opts):
    "Args: list of .ufl filenames."
    pass #_add_skipping_options(opts)

def run_str(options, args):
    "Print str of the objects found in .ufl file(s)."
    load_and_print_objects(args, str,
                           skip_elements=options.skip_elements,
                           skip_coefficients=options.skip_coefficients,
                           skip_expressions=options.skip_expressions,
                           skip_forms=options.skip_forms)
    return 0


def add_repr_options(opts):
    "Args: list of .ufl filenames."
    pass #_add_skipping_options(opts)

def run_repr(options, args):
    "Print repr of the objects found in .ufl file(s)."
    load_and_print_objects(args, repr,
                           skip_elements=options.skip_elements,
                           skip_coefficients=options.skip_coefficients,
                           skip_expressions=options.skip_expressions,
                           skip_forms=options.skip_forms)
    return 0


def add_tree_options(opts):
    "Args: list of .ufl filenames."
    pass
    #opts.add_option('--skip-coefficients', action='store_true', default=False)
    #opts.add_option('--skip-expressions', action='store_true', default=False)
    #opts.add_option('--skip-forms', action='store_true', default=False)

def run_tree(options, args):
    "Print tree_format of the objects found in .ufl file(s)."
    from ufl.algorithms import tree_format
    load_and_print_objects(args, tree_format,
                           skip_elements=True,
                           skip_coefficients=options.skip_coefficients,
                           skip_expressions=options.skip_expressions,
                           skip_forms=options.skip_forms)
    return 0
