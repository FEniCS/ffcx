# Copyright (C) 2009-2017 Anders Logg, Martin Sandve Alnæs, Marie E. Rognes,
# Kristian B. Oelgaard, and others
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Compiler stage 4: Code generation.

This module implements the generation of C code for the body of each
UFCx function from an intermediate representation (IR).
"""

from __future__ import annotations

import logging
import typing
from importlib import import_module

import numpy.typing as npt
import ufl

from ffcx.ir.representation import DataIR
from ffcx.options import get_language

logger = logging.getLogger("ffcx")


class CodeBlocks(typing.NamedTuple):
    """Storage of code blocks of the form (declaration, implementation).

    Blocks for the functions of external operators, integrals, forms and
    expressions, and start and end of file output. The functions of external
    operators come first, ahead of the kernels that call them.
    """

    file_pre: list[tuple[str, str]]
    external_functions: list[tuple[str, str]]
    expressions: list[tuple[str, str]]
    integrals: list[tuple[str, str]]
    forms: list[tuple[str, str]]
    file_post: list[tuple[str, str]]


def generate_code(
    ir: DataIR, options: dict[str, int | float | npt.DTypeLike]
) -> tuple[CodeBlocks, tuple[str, ...]]:
    """Generate code blocks from intermediate representation."""
    logger.info(79 * "*")
    logger.info("Compiler stage 3: Generating code")
    logger.info(79 * "*")

    mod = import_module(get_language(options))

    integral_generator = mod.integral.generator
    form_generator = mod.form.generator
    expression_generator = mod.expression.generator
    file_generator = mod.file.generator

    # First, as the backend checks that it can call the functions of the
    # external operators before the kernels call them
    externals = _external_operators(ir)
    if externals and not hasattr(mod, "external"):
        raise NotImplementedError(
            f"External operators are not supported by the {options['language']} backend."
        )
    code_external = mod.external.generator(externals, options) if externals else []
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
        external_functions=code_external,
        integrals=code_integrals,
        forms=code_forms,
        expressions=code_expressions,
        file_post=[code_file_post],
    ), mod.file.suffixes


def _external_operators(ir: DataIR) -> list[ufl.ExternalOperator]:
    """Get the external operators the kernels evaluate, each once."""
    sub_kernels = [k for i in ir.integrals for k in i.sub_kernels]
    sub_kernels += [k for e in ir.expressions for k in e.sub_kernels]
    externals = (getattr(k.target, "external", None) for k in sub_kernels)
    return list(dict.fromkeys(e for e in externals if e is not None))
