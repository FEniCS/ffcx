
from uflacs.analysis.graph import build_graph

def testelementmapping(form): # TODO: Maybe modify and add as ufl unit test?
    from ufl import quadrilateral
    from ufl.algorithms import extract_elements

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

def test_build_graph(form_data):
    form = form_data.preprocessed_form
    for integrand in [itg.integrand() for itg in form.integrals()]:
        G = build_graph([integrand], DEBUG=1)
        print G.total_value_size, G.total_unique_symbols

# TODO: Maybe move this to a utility file and reuse. Maybe make iterator versions over forms/expressions?
def for_each_ufl_file(filenames, file_cb=None, form_cb=None, form_data_cb=None):
    from ufl.algorithms import load_ufl_file
    for fn in filenames:
        data = load_ufl_file(fn)
        if file_cb:
            file_cb(data)
        for form in data.forms:
            if form_cb:
                form_cb(form)
            if form_data_cb:
                form_data_cb(form.compute_form_data())

def add_debug_options(opts):
    "Args: list of .ufl filenames."
    pass

def run_debug(options, args):
    "For use while developing. Run whatever debugging code has currently been implemented here."
    #for_each_ufl_file(args, form_cb=testelementmapping)
    for_each_ufl_file(args, form_data_cb=test_build_graph)
