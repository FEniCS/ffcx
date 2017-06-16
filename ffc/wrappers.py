# -*- coding: utf-8 -*-

# Copyright (C) 2010-2017 Anders Logg
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

# Python modules
from itertools import chain

# FFC modules
from ffc.log import begin, end, info, error
from ffc.utils import all_equal
from ffc.classname import make_classname
from ffc.backends.dolfin.wrappers import generate_dolfin_code
from ffc.backends.dolfin.capsules import UFCElementNames, UFCFormNames

__all__ = ["generate_wrapper_code"]

# FIXME: More clean-ups needed here.


def generate_wrapper_code(analysis, prefix, object_names, parameters):
    "Generate code for additional wrappers."

    # Skip if wrappers not requested
    if not parameters["format"] == "dolfin":
        return None

    # Return dolfin wrapper
    return _generate_dolfin_wrapper(analysis, prefix, object_names, parameters)


def _generate_dolfin_wrapper(analysis, prefix, object_names, parameters):

    begin("Compiler stage 4.1: Generating additional wrapper code")

    # Encapsulate data
    (capsules, common_space) = _encapsulate(prefix, object_names, analysis, parameters)

    # Generate code
    info("Generating wrapper code for DOLFIN")
    code = generate_dolfin_code(prefix, "",
                                capsules, common_space,
                                error_control=parameters["error_control"])
    code += "\n\n"
    end()

    return code


def _encapsulate(prefix, object_names, analysis, parameters):

    # Extract data from analysis
    form_datas, elements, element_map, domains = analysis

    # FIXME: Encapsulate domains?

    num_form_datas = len(form_datas)
    common_space = False

    # Special case: single element
    if num_form_datas == 0:
        capsules = _encapsule_element(prefix, elements)

    # Special case: with error control
    elif parameters["error_control"] and num_form_datas == 11:
        capsules = [_encapsule_form(prefix, object_names, form_data, i, element_map)
                    for (i, form_data) in enumerate(form_datas[:num_form_datas - 1])]
        capsules += [_encapsule_form(prefix, object_names, form_datas[-1], num_form_datas - 1,
                                     element_map, "GoalFunctional")]
    # Otherwise: generate standard capsules for each form
    else:
        capsules = [_encapsule_form(prefix, object_names, form_data, i, element_map) for
                    (i, form_data) in enumerate(form_datas)]
        # Check if all argument elements are equal
        elements = []
        for form_data in form_datas:
            elements += form_data.argument_elements
        common_space = all_equal(elements)

    return (capsules, common_space)


def _encapsule_form(prefix, object_names, form_data, i, element_map, superclassname=None):
    element_numbers = [element_map[e] for e in form_data.argument_elements + form_data.coefficient_elements]

    if superclassname is None:
        superclassname = "Form"

    form_names = UFCFormNames(
        object_names.get(id(form_data.original_form), "%d" % i),
        [object_names.get(id(obj), "w%d" % j) for j, obj in enumerate(form_data.reduced_coefficients)],
        make_classname(prefix, "form", i),
        [make_classname(prefix, "finite_element", j) for j in element_numbers],
        [make_classname(prefix, "dofmap", j) for j in element_numbers],
        superclassname)

    return form_names


def _encapsule_element(prefix, elements):
    element_number = len(elements) - 1  # eh? this doesn't make any sense
    args = ("0",
            [make_classname(prefix, "finite_element", element_number)],
            [make_classname(prefix, "dofmap", element_number)])
    return UFCElementNames(*args)
