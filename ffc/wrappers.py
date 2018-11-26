# -*- coding: utf-8 -*-
# Copyright (C) 2010-2017 Anders Logg
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import logging

from ffc.backends import dolfin

logger = logging.getLogger(__name__)


def generate_wrapper_code(analysis, prefix, object_names, classnames, parameters):

    logger.info("Compiler stage 4.1: Generating additional wrapper code")

    # Encapsulate data
    capsules, common_space = _encapsulate(prefix, object_names, classnames, analysis, parameters)

    # Generate code
    logger.info("Generating wrapper code for DOLFIN")
    code_h, code_c = dolfin.generate_wrappers(prefix, capsules, common_space)
    code_h += "\n"
    code_c += "\n"

    return code_h, code_c


def _encapsulate(prefix, object_names, classnames, analysis, parameters):

    # Extract data from analysis
    form_data, elements, element_map, domains = analysis

    # FIXME: Encapsulate domains?

    if not form_data:
        return _encapsulate_elements(elements, object_names, classnames), False
    else:
        return _encapsule_forms(prefix, object_names, classnames, form_data, element_map)


def _encapsulate_elements(elements, object_names, classnames):
    """Generate capsules for each element named in the input (no wrappers
    for subelements will be created)"""
    assert not classnames["coordinate_maps"], "Need to fix element wrappers for coordinate mappings."

    capsules = []
    for i, element in enumerate(elements):
        name = object_names.get(id(element), None)
        if name:
            args = {
                "name": name,
                "element_classname": classnames["elements"][i],
                "dofmap_classname": classnames["dofmaps"][i],
                "coordinate_mapping_classname": None
            }
            capsules.append(dolfin.UFCElementNames(**args))

    return capsules


def _encapsule_forms(prefix, object_names, classnames, form_data, element_map):

    # FIXME: Figure what to do with coordinate maps. Can there be more
    # than 1?
    assert len(classnames["coordinate_maps"]) == 1

    capsules = []
    for i, form in enumerate(form_data):
        element_numbers = [
            element_map[e] for e in form.argument_elements + form.coefficient_elements
        ]

        name = object_names.get(id(form.original_form), "%d" % i)
        coefficient_names = [
            object_names.get(id(obj), "w%d" % j) for j, obj in enumerate(form.reduced_coefficients)
        ]
        ufc_form_name = classnames["forms"][i]
        ufc_elements = [classnames["elements"][j] for j in element_numbers]
        ufc_dofmaps = [classnames["dofmaps"][j] for j in element_numbers]
        ufc_cmaps = [classnames["coordinate_maps"][0]]

        capsules.append(
            dolfin.UFCFormNames(name, coefficient_names, ufc_form_name, ufc_elements, ufc_dofmaps,
                                ufc_cmaps))

    # Build list of all argument elements to which if all are equal
    elements = [element for form in form_data for element in form.argument_elements]

    return capsules, all(x == elements[0] for x in elements)
