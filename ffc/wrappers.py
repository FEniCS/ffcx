__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2010-01-18"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2011-01-26

# Python modules
from itertools import chain

# FFC modules
from ffc.log import begin, end, info, error
from ffc.utils import all_equal
from ffc.cpp import format

__all__ = ["generate_wrapper_code"]

# FIXME: More clean-ups needed here.

def generate_wrapper_code(analysis, prefix, parameters):
    "Generate code for additional wrappers."

    # Skip if wrappers not requested
    if not parameters["format"] == "dolfin":
        return None

    # Check that we can import wrappers from dolfin
    try:
        import dolfin_utils.wrappers
    except:
        error("Unable to generate new DOLFIN wrappers, missing module dolfin_utils.wrappers.")

    # Return dolfin wrapper
    return _generate_dolfin_wrapper(analysis, prefix, parameters)

def _generate_dolfin_wrapper(analysis, prefix, parameters):

    begin("Compiler stage 4.1: Generating additional wrapper code")

    # Extract data from analysis
    forms, elements, element_map = analysis

    # Encapsulate data
    (capsules, common_space) = _encapsulate(prefix, forms, elements,
                                           element_map, parameters)

    # Generate code
    info("Generating wrapper code for DOLFIN")
    from dolfin_utils.wrappers import generate_dolfin_code
    code = generate_dolfin_code(prefix, "", capsules, common_space) + "\n\n"
    end()

    return code

def _encapsulate(prefix, forms, elements, element_map, parameters):

    num_forms = len(forms)
    common_space = False

    # Special case: single element
    if num_forms == 0:
        capsules = encapsule_element(prefix, elements)

    # Special case: with error control
    elif (parameters["error_control"] and num_forms == 11):
        capsules = [_encapsule_form(prefix, form, i, element_map) for
                    (i, form) in enumerate(forms[:num_forms-1])]
        capsules += [_encapsule_form(prefix, forms[-1], num_forms-1,
                                     element_map, "GoalFunctional")]

    # Otherwise: generate standard capsules for each form
    else:
        capsules = [_encapsule_form(prefix, form, i, element_map) for
                    (i, form) in enumerate(forms)]

        # Check if all elements are equal
        elements = []
        for form in forms:
            elements += form.form_data().elements[:form.form_data().rank]
        common_space = all_equal(elements)

    return (capsules, common_space)


def _encapsule_form(prefix, form, i, element_map, superclassname=None):
    element_numbers = [element_map[e] for e in form.form_data().elements]

    if superclassname is None:
        superclassname = "Form"

    from dolfin_utils.wrappers import UFCFormNames
    form_names = UFCFormNames("%d" % i,
                              form.form_data().coefficient_names,
                              format["classname form"](prefix, i),
                              [format["classname finite_element"](prefix, j)
                               for j in element_numbers],
                              [format["classname dof_map"](prefix, j)
                               for j in element_numbers],
                              superclassname)
    return form_names

def _encapsule_element(prefix, elements):
    element_number = len(elements) - 1
    args = ("0",
            [format["classname finite_element"](prefix, element_number)],
            [format["classname dof_map"](prefix, element_number)])
    from dolfin_utils.wrappers import UFCElementNames
    return UFCElementNames(*args)
