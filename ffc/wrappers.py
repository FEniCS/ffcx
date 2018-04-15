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
from ffc.backends.dolfin import generate_dolfin_code
from ffc.backends.dolfin.capsules import UFCElementNames, UFCFormNames

__all__ = ["generate_wrapper_code"]

# FIXME: More clean-ups needed here.


def generate_wrapper_code(analysis, prefix, object_names, classnames,
                          parameters):
    "Generate code for additional wrappers."

    # Skip if wrappers not requested
    if not parameters["format"] == "dolfin":
        return None

    # Return dolfin wrapper
    return _generate_dolfin_wrapper(analysis, prefix, object_names, classnames,
                                    parameters)


def _generate_dolfin_wrapper(analysis, prefix, object_names, classnames,
                             parameters):

    begin("Compiler stage 4.1: Generating additional wrapper code")

    # Encapsulate data
    capsules, common_space = _encapsulate(prefix, object_names, classnames,
                                          analysis, parameters)

    # Generate code
    info("Generating wrapper code for DOLFIN")
    code_h, code_c = generate_dolfin_code(prefix, capsules, common_space)
    code_h += "\n\n"
    code_c += "\n\n"
    end()

    return code_h, code_c


def _encapsulate(prefix, object_names, classnames, analysis, parameters):

    # Extract data from analysis
    form_datas, elements, element_map, domains = analysis

    # FIXME: Encapsulate domains?

    if not form_datas:
        # Generate capsules for each element
        common_space = False
        capsules = []
        for i in range(len(elements)):
            assert len(classnames["coordinate_maps"]) == 0
            element_number = i
            args = {
                "name": i,
                "element_classname": classnames["elements"][i],
                "dofmap_classname": classnames["dofmaps"][i],
                "coordinate_mapping_classname": None
            }
            capsules.append(UFCElementNames(**args))
    else:
        # Generate capsules for each form
        capsules = [
            _encapsule_form(prefix, object_names, classnames, form_data, i,
                            element_map)
            for (i, form_data) in enumerate(form_datas)
        ]
        # Check if all argument elements are equal
        elements = []
        for form_data in form_datas:
            elements += form_data.argument_elements
        common_space = all_equal(elements)

    return capsules, common_space


def _encapsule_form(prefix, object_names, classnames, form_data, i,
                    element_map):
    element_numbers = [
        element_map[e]
        for e in form_data.argument_elements + form_data.coefficient_elements
    ]

    # FIXME: Figure what to do with coordinate maps. Can there be more than 1?
    assert len(classnames["coordinate_maps"]) == 1
    form_names = UFCFormNames(
        object_names.get(id(form_data.original_form), "%d" % i), [
            object_names.get(id(obj), "w%d" % j)
            for j, obj in enumerate(form_data.reduced_coefficients)
        ], classnames["forms"][i],
        [classnames["elements"][j] for j in element_numbers],
        [classnames["dofmaps"][j]
         for j in element_numbers], [classnames["coordinate_maps"][0]])

    return form_names
