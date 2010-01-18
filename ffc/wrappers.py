__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2010-01-18"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-01-18

# Python modules
from itertools import chain

# FFC modules
from ffc.log import begin, end, info, error
from ffc.utils import all_equal
from ffc.cpp import format

def generate_wrapper_code(analysis, prefix, options):
    "Generate code for additional wrappers."

    # Skip if wrappers not requested
    if not options["format"] == "dolfin": return None

    # Try importing DOLFIN wrapper utils
    try:
        from dolfin_utils.wrappers import generate_dolfin_code, UFCFormNames
    except:
        error("Unable to generate DOLFIN wrappers, missing module dolfin_utils.wrappers.")

    begin("Compiler stage 4.1: Generating additional wrapper code")

    # Extract data from analysis
    form_and_data, elements, element_map = analysis

    # Special case: single element
    #if len(generated_forms) == 1 and generated_forms[0][1].form is None:
    #    fn = UFCFormNames("0",
    #                      [],
    #                      format["classname form"](prefix, 0),
    #                      [format["classname finite_element"](prefix, 0, (0,))],
    #                      [format["classname dof_map"](prefix, 0, (0,))])
    #    return generate_dolfin_code(prefix, "", [fn], (0, 0), False) + "\n\n"

    # Generate name data for each form
    form_names = []
    for (i, (form, form_data)) in enumerate(form_and_data):
        element_numbers = [element_map[e] for e in form_data.elements]

        print form_data.coefficient_names

        form_names.append(UFCFormNames("%d" % i,
                                       form_data.coefficient_names,
                                       format["classname form"](prefix, i),
                                       [format["classname finite_element"](prefix, j) for j in element_numbers],
                                       [format["classname dof_map"](prefix, j) for j in element_numbers]))

    # Check if all elements are equal
    if all_equal(elements):
        common_space = (0, 0)
    else:
        common_space = None

    # Generate code
    info("Generating wrapper code for DOLFIN")
    code = generate_dolfin_code(prefix, "", form_names, common_space, False) + "\n\n"

    end()

    return code
