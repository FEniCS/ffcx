
from uflacs.analysis.graph import build_graph

def test_build_graph(form_data):
    form = form_data.preprocessed_form
    for integrand in [itg.integrand() for itg in form.integrals()]:
        G = build_graph([integrand], DEBUG=1)
        print(G.total_value_size, G.total_unique_symbols)

# TODO: Maybe move this to a utility file and reuse. Maybe make iterator versions over forms/expressions?
def for_each_ufl_file(filenames, file_cb=None, form_cb=None, form_data_cb=None):
    from ufl.algorithms import load_ufl_file, compute_form_data
    for fn in filenames:
        data = load_ufl_file(fn)
        if file_cb:
            file_cb(data)
        for form in data.forms:
            if form_cb:
                form_cb(form)
            if form_data_cb:
                form_data_cb(compute_form_data(form))

def add_debug_options(opts):
    "Args: list of .ufl filenames."
    pass

def run_debug(options, args):
    "For use while developing. Run whatever debugging code has currently been implemented here."
    for_each_ufl_file(args, form_data_cb=test_build_graph)
