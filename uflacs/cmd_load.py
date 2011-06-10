
from uflacs.utils.log import info

def add_load_options(opts):
    "Args: list of .ufl filenames."
    pass

def run_load(options, args):
    "Load .ufl file(s) and print the members of the loaded data structure."
    from uflacs.utils.str_utils import format_members
    from ufl.algorithms import load_ufl_file
    for fn in args:
        info("Loading file '%s'..." % (fn,))
        data = load_ufl_file(fn)
        info(format_members(data))
    return 0
