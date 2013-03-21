# Copyright (C) 2010-2013 Anders Logg
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.
#
# First added:  2010-01-18
# Last changed: 2013-03-21

# Python modules
from itertools import chain

# FFC modules
from ffc.log import begin, end, info, error
from ffc.utils import all_equal
from ffc.cpp import format
from ffc.dolfin.wrappers import generate_dolfin_code
from ffc.dolfin.capsules import UFCElementNames, UFCFormNames

__all__ = ["generate_wrapper_code"]

# FIXME: More clean-ups needed here.

def generate_wrapper_code(analysis, prefix, parameters):
    "Generate code for additional wrappers."

    # Skip if wrappers not requested
    if not parameters["format"] == "dolfin":
        return None

    # Return dolfin wrapper
    return _generate_dolfin_wrapper(analysis, prefix, parameters)

def _generate_dolfin_wrapper(analysis, prefix, parameters):

    begin("Compiler stage 4.1: Generating additional wrapper code")

    # Encapsulate data
    (capsules, common_space) = _encapsulate(prefix, analysis, parameters)

    # Generate code
    info("Generating wrapper code for DOLFIN")
    code = generate_dolfin_code(prefix, "",
                                capsules, common_space,
                                error_control=parameters["error_control"])
    code += "\n\n"
    end()

    return code

def _encapsulate(prefix, analysis, parameters):

    # Extract data from analysis
    form_datas, elements, element_map = analysis

    num_form_datas = len(form_datas)
    common_space = False

    # Special case: single element
    if num_form_datas == 0:
        capsules = _encapsule_element(prefix, elements)

    # Special case: with error control
    elif (parameters["error_control"] and num_form_datas == 11):
        capsules = [_encapsule_form(prefix, form_data, i, element_map) for
                    (i, form_data) in enumerate(form_datas[:num_form_datas-1])]
        capsules += [_encapsule_form(prefix, form_datas[-1], num_form_datas-1,
                                     element_map, "GoalFunctional")]

    # Otherwise: generate standard capsules for each form
    else:
        capsules = [_encapsule_form(prefix, form_data, i, element_map) for
                    (i, form_data) in enumerate(form_datas)]

        # Check if all elements are equal
        elements = []
        for form_data in form_datas:
            elements += form_data.elements[:form_data.rank]
        common_space = all_equal(elements)

    return (capsules, common_space)


def _encapsule_form(prefix, form_data, i, element_map, superclassname=None):
    element_numbers = [element_map[e] for e in form_data.elements]

    if superclassname is None:
        superclassname = "Form"

    form_names = UFCFormNames(form_data.name or "%d" % i,
                              form_data.coefficient_names,
                              format["classname form"](prefix, i),
                              [format["classname finite_element"](prefix, j)
                               for j in element_numbers],
                              [format["classname dofmap"](prefix, j)
                               for j in element_numbers],
                              superclassname)
    return form_names

def _encapsule_element(prefix, elements):
    element_number = len(elements) - 1
    args = ("0",
            [format["classname finite_element"](prefix, element_number)],
            [format["classname dofmap"](prefix, element_number)])
    return UFCElementNames(*args)
