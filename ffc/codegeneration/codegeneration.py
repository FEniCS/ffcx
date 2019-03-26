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

import logging
from collections import namedtuple

from ffc.codegeneration.coordinate_mapping import \
    ufc_coordinate_mapping_generator
from ffc.codegeneration.dofmap import ufc_dofmap_generator
from ffc.codegeneration.finite_element import \
    generator as ufc_finite_element_generator
from ffc.codegeneration.form import ufc_form_generator
from ffc.codegeneration.integrals import ufc_integral_generator

logger = logging.getLogger(__name__)

code_blocks = namedtuple('code_blocks', ['elements', 'dofmaps',
                                         'coordinate_mappings', 'integrals', 'forms'])


def generate_code(analysis, object_names, ir, parameters):
    """Generate code from intermediate representation."""

    logger.debug("Compiler stage 4: Generating code")

    # Generate code for finite_elements
    logger.debug("Generating code for {} finite_element(s)".format(len(ir.elements)))
    code_finite_elements = [
        ufc_finite_element_generator(ir, parameters) for ir in ir.elements
    ]

    # Generate code for dofmaps
    logger.debug("Generating code for {} dofmap(s)".format(len(ir.dofmaps)))
    code_dofmaps = [ufc_dofmap_generator(ir, parameters) for ir in ir.dofmaps]

    # Generate code for coordinate_mappings
    logger.debug("Generating code for {} coordinate_mapping(s)".format(len(ir.coordinate_mappings)))
    code_coordinate_mappings = [
        ufc_coordinate_mapping_generator(ir, parameters) for ir in ir.coordinate_mappings
    ]

    # Generate code for integrals
    logger.debug("Generating code for integrals")
    code_integrals = [ufc_integral_generator(ir, parameters) for ir in ir.integrals]

    # Generate code for forms
    logger.debug("Generating code for forms")
    # FIXME: add coefficient names to IR
    coefficient_names = []
    for form in analysis.form_data:
        names = [
            object_names.get(id(obj), "w%d" % j) for j, obj in enumerate(form.reduced_coefficients)
        ]
        coefficient_names.append(names)
    code_forms = [ufc_form_generator(ir, cnames, parameters) for ir, cnames in zip(ir.forms, coefficient_names)]

    return code_blocks(elements=code_finite_elements, dofmaps=code_dofmaps,
                       coordinate_mappings=code_coordinate_mappings, integrals=code_integrals,
                       forms=code_forms)
