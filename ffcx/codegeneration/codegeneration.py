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
import typing
import sys
from importlib import import_module

logger = logging.getLogger("ffcx")


class CodeBlocks(typing.NamedTuple):
    """
    Storage of code blocks of the form (declaration, implementation).

    Blocks for elements, dofmaps, integrals, forms and expressions,
    and start and end of file output
    """

    file_pre: typing.List[typing.Tuple[str, str]]
    elements: typing.List[typing.Tuple[str, str]]
    dofmaps: typing.List[typing.Tuple[str, str]]
    integrals: typing.List[typing.Tuple[str, str]]
    forms: typing.List[typing.Tuple[str, str]]
    expressions: typing.List[typing.Tuple[str, str]]
    file_post: typing.List[typing.Tuple[str, str]]


def generate_code(ir, options) -> typing.Tuple[CodeBlocks, typing.Tuple[str]]:
    """Generate code blocks from intermediate representation."""
    logger.info(79 * "*")
    logger.info("Compiler stage 3: Generating code")
    logger.info(79 * "*")

    lang = options.get("language", "C")
    try:
        # Built-in
        mod = import_module(f"ffcx.codegeneration.{lang}")
    except ImportError:
        # User defined, in current directory
        sys.path.append(".")
        mod = import_module(f"{lang}")

    finite_element_generator = mod.finite_element.generator
    dofmap_generator = mod.dofmap.generator
    integral_generator = mod.integrals.generator
    form_generator = mod.form.generator
    expression_generator = mod.expressions.generator
    file_generator = mod.file.generator

    # Generate code for finite_elements
    code_finite_elements = [finite_element_generator(element_ir, options) for element_ir in ir.elements]
    code_dofmaps = [dofmap_generator(dofmap_ir, options) for dofmap_ir in ir.dofmaps]
    code_integrals = [integral_generator(integral_ir, options) for integral_ir in ir.integrals]
    code_forms = [form_generator(form_ir, options) for form_ir in ir.forms]
    code_expressions = [expression_generator(expression_ir, options) for expression_ir in ir.expressions]
    code_file_pre, code_file_post = file_generator(options)
    return CodeBlocks(file_pre=[code_file_pre], elements=code_finite_elements, dofmaps=code_dofmaps,
                      integrals=code_integrals, forms=code_forms, expressions=code_expressions,
                      file_post=[code_file_post]), mod.suffixes
