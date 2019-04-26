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

from ffc.codegeneration.finite_element import generator as finite_element_generator
from ffc.codegeneration.coordinate_mapping import \
    generator as coordinate_mapping_generator
from ffc.codegeneration.dofmap import generator as dofmap_generator
from ffc.codegeneration.form import generator as form_generator
from ffc.codegeneration.integrals import generator as integral_generator

logger = logging.getLogger(__name__)

code_blocks = namedtuple('code_blocks', ['elements', 'dofmaps',
                                         'coordinate_mappings', 'integrals', 'forms'])


def generate_code(ir, parameters):
    """Generate code blocks from intermediate representation.

    """

    logger.debug("Compiler stage 4: Generating code")

    # Generate code for finite_elements
    logger.debug("Generating code for {} finite_element(s)".format(len(ir.elements)))
    code_finite_elements = [finite_element_generator(element_ir, parameters) for element_ir in ir.elements]

    # Generate code for dofmaps
    logger.debug("Generating code for {} dofmap(s)".format(len(ir.dofmaps)))
    code_dofmaps = [dofmap_generator(dofmap_ir, parameters) for dofmap_ir in ir.dofmaps]

    # Generate code for coordinate_mappings
    logger.debug("Generating code for {} coordinate_mapping(s)".format(len(ir.coordinate_mappings)))
    code_coordinate_mappings = [coordinate_mapping_generator(cmap_ir, parameters) for cmap_ir in ir.coordinate_mappings]

    # Generate code for integrals
    logger.debug("Generating code for integrals")
    code_integrals = [integral_generator(integral_ir, parameters) for integral_ir in ir.integrals]

    # Generate code for forms
    logger.debug("Generating code for forms")
    code_forms = [form_generator(form_ir, parameters) for form_ir in ir.forms]

    return code_blocks(elements=code_finite_elements, dofmaps=code_dofmaps,
                       coordinate_mappings=code_coordinate_mappings, integrals=code_integrals,
                       forms=code_forms)
