# Copyright (C) 2009-2017 Anders Logg, Martin Sandve Alnæs, Marie E. Rognes,
# Kristian B. Oelgaard, and others
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Compiler stage 4: Code generation.

This module implements the generation of C code for the body of each
UFC function from an intermediate representation (IR).

"""

import logging
from collections import namedtuple

from ffcx.codegeneration.dofmap import generator as dofmap_generator
from ffcx.codegeneration.expressions import generator as expression_generator
from ffcx.codegeneration.finite_element import \
    generator as finite_element_generator
from ffcx.codegeneration.form import generator as form_generator
from ffcx.codegeneration.integrals import generator as integral_generator

logger = logging.getLogger("ffcx")

code_blocks = namedtuple("code_blocks", ["elements", "dofmaps", "integrals", "forms", "expressions"])


def generate_code(ir, parameters):
    """Generate code blocks from intermediate representation."""
    logger.info(79 * "*")
    logger.info("Compiler stage 3: Generating code")
    logger.info(79 * "*")

    # Generate code for finite_elements
    code_finite_elements = [finite_element_generator(element_ir, parameters) for element_ir in ir.elements]
    code_dofmaps = [dofmap_generator(dofmap_ir, parameters) for dofmap_ir in ir.dofmaps]
    code_integrals = [integral_generator(integral_ir, parameters) for integral_ir in ir.integrals]
    code_forms = [form_generator(form_ir, parameters) for form_ir in ir.forms]
    code_expressions = [expression_generator(expression_ir, parameters) for expression_ir in ir.expressions]

    return code_blocks(elements=code_finite_elements, dofmaps=code_dofmaps,
                       integrals=code_integrals, forms=code_forms, expressions=code_expressions)
