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

from __future__ import annotations

import logging
import sys
import typing
from importlib import import_module

import numpy.typing as npt

from ffcx.ir.representation import DataIR

logger = logging.getLogger("ffcx")


class CodeBlocks(typing.NamedTuple):
    """Storage of code blocks of the form (declaration, implementation).

    Blocks for integrals, forms and expressions, and start and end of file output
    """

    file_pre: list[tuple[str, str]]
    integrals: list[tuple[str, str]]
    forms: list[tuple[str, str]]
    expressions: list[tuple[str, str]]
    file_post: list[tuple[str, str]]


def generate_code(
    ir: DataIR, options: dict[str, int | float | npt.DTypeLike]
) -> tuple[CodeBlocks, tuple[str, str]]:
    """Generate code blocks from intermediate representation."""
    logger.info(79 * "*")
    logger.info("Compiler stage 3: Generating code")
    logger.info(79 * "*")

    lang = options.get("language", "C")
    try:
        # Built-in
        mod = import_module(f"ffcx.codegeneration.{lang}")
    except ImportError:
        # User defined language (experimental)
        store_path = sys.path
        sys.path = ["."]
        mod = import_module(f"{lang}")
        sys.path = store_path

    integral_generator = mod.integrals.generator
    form_generator = mod.form.generator
    expression_generator = mod.expressions.generator
    file_generator = mod.file.generator

    code_integrals = [
        integral_generator(integral_ir, domain, options)
        for integral_ir in ir.integrals
        for domain in set(i[0] for i in integral_ir.expression.integrand.keys())
    ]
    code_forms = [form_generator(form_ir, options) for form_ir in ir.forms]
    code_expressions = [
        expression_generator(expression_ir, options) for expression_ir in ir.expressions
    ]
    code_file_pre, code_file_post = file_generator(options)
    return CodeBlocks(
        file_pre=[code_file_pre],
        integrals=code_integrals,
        forms=code_forms,
        expressions=code_expressions,
        file_post=[code_file_post],
    ), mod.suffixes
