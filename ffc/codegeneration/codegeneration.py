# -*- coding: utf-8 -*-
# Copyright (C) 2009-2017 Anders Logg, Martin Sandve Alnæs, Marie E. Rognes,
# Kristian B. Oelgaard, and others
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Compiler stage 4: Code generation

This module implements the generation of C code for the body of each
UFC function from an intermediate representation (IR).

"""

import itertools
import logging

from ffc.codegeneration.coordinate_mapping import ufc_coordinate_mapping_generator
from ffc.codegeneration.dofmap import ufc_dofmap_generator
from ffc.codegeneration.finite_element import generator as ufc_finite_element_generator
from ffc.codegeneration.form import ufc_form_generator
from ffc.codegeneration.integrals import ufc_integral_generator

logger = logging.getLogger(__name__)


def generate_code(analysis, object_names, ir, parameters, jit):
    """Generate code from intermediate representation."""

    logger.debug("Compiler stage 4: Generating code")

    full_ir = ir

    # FIXME: This has global side effects
    # Set code generation parameters
    # set_float_formatting(parameters["precision"])
    # set_exception_handling(parameters["convert_exceptions_to_warnings"])

    # Extract data from analysis
    form_data, elements, element_map, domains = analysis

    # Extract representations
    ir_finite_elements, ir_dofmaps, ir_coordinate_mappings, ir_integrals, ir_forms = ir

    # Generate code for finite_elements
    logger.debug("Generating code for {} finite_element(s)".format(len(ir_finite_elements)))
    code_finite_elements = [
        ufc_finite_element_generator(ir, parameters) for ir in ir_finite_elements
    ]

    # Generate code for dofmaps
    logger.debug("Generating code for {} dofmap(s)".format(len(ir_dofmaps)))
    code_dofmaps = [ufc_dofmap_generator(ir, parameters) for ir in ir_dofmaps]

    # Generate code for coordinate_mappings
    logger.debug("Generating code for {} coordinate_mapping(s)".format(len(ir_coordinate_mappings)))
    code_coordinate_mappings = [
        ufc_coordinate_mapping_generator(ir, parameters) for ir in ir_coordinate_mappings
    ]

    # Generate code for integrals
    logger.debug("Generating code for integrals")
    code_integrals = [ufc_integral_generator(ir, parameters) for ir in ir_integrals]

    # Generate code for forms
    logger.debug("Generating code for forms")
    # FIXME: add coefficient names it IR
    coefficient_names = []
    for form in form_data:
        names = [
            object_names.get(id(obj), "w%d" % j) for j, obj in enumerate(form.reduced_coefficients)
        ]
        coefficient_names.append(names)
    code_forms = [ufc_form_generator(ir, cnames, parameters) for ir, cnames in zip(ir_forms, coefficient_names)]

    # Extract additional includes
    includes = _extract_includes(full_ir, code_integrals, jit)

    return (code_finite_elements, code_dofmaps, code_coordinate_mappings, code_integrals,
            code_forms, includes)


def _extract_includes(full_ir, code_integrals, jit):
    ir_finite_elements, ir_dofmaps, ir_coordinate_mappings, ir_integrals, ir_forms = full_ir

    # Includes added by representations
    includes = set()
    # for code in code_integrals:
    #     includes.update(code["additional_includes_set"])

    # Includes for dependencies in jit mode
    if jit:
        dep_includes = set()
        for ir in ir_finite_elements:
            dep_includes.update(_finite_element_jit_includes(ir))
        for ir in ir_dofmaps:
            dep_includes.update(_dofmap_jit_includes(ir))
        for ir in ir_coordinate_mappings:
            dep_includes.update(_coordinate_mapping_jit_includes(ir))
        # for ir in ir_integrals:
        #    dep_includes.update(_integral_jit_includes(ir))
        for ir in ir_forms:
            dep_includes.update(_form_jit_includes(ir))
        includes.update(['#include "{}"'.format(inc) for inc in dep_includes])

    return includes


def _finite_element_jit_includes(ir):
    classnames = ir["create_sub_element"]
    postfix = "_finite_element"
    return [classname.rpartition(postfix)[0] + ".h" for classname in classnames]


def _dofmap_jit_includes(ir):
    classnames = ir["create_sub_dofmap"]
    postfix = "_dofmap"
    return [classname.rpartition(postfix)[0] + ".h" for classname in classnames]


def _coordinate_mapping_jit_includes(ir):
    classnames = [
        ir["coordinate_finite_element_classname"], ir["scalar_coordinate_finite_element_classname"]
    ]
    postfix = "_finite_element"
    return [classname.rpartition(postfix)[0] + ".h" for classname in classnames]


def _form_jit_includes(ir):
    # Gather all header names for classes that are separately compiled
    # For finite_element and dofmap the module and header name is the prefix,
    # extracted here with .split, and equal for both classes so we skip dofmap here:
    classnames = list(
        itertools.chain(ir["create_finite_element"], ir["create_coordinate_finite_element"]))
    postfix = "_finite_element"
    includes = [classname.rpartition(postfix)[0] + ".h" for classname in classnames]

    classnames = ir["create_coordinate_mapping"]
    postfix = "_coordinate_mapping"
    includes += [classname.rpartition(postfix)[0] + ".h" for classname in classnames]
    return includes
