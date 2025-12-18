# Copyright (C) 2009-2018 Anders Logg and Garth N. Wells
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Compiler stage 5: Code formatting.

This module implements the formatting of UFC code from a given
dictionary of generated C++ code for the body of each UFC function.

It relies on templates for UFC code available as part of the module
ufcx_utils.

"""

from __future__ import annotations

import logging
from pathlib import Path

from ffcx.codegeneration.codegeneration import CodeBlocks

logger = logging.getLogger("ffcx")


def format_code(code_blocks: CodeBlocks) -> list[str]:
    """Format given code in UFC format. Returns two strings with header and source file contents."""
    logger.info(79 * "*")
    logger.info("Compiler stage 5: Formatting code")
    logger.info(79 * "*")

    code = [""] * len(code_blocks[0][0])

    for block in code_blocks:
        for i in range(len(code)):
            code[i] += "".join([c[i] for c in block])

    return code


def write_code(code: list[str], prefix: str, suffixes: tuple[str, ...], output_dir: str) -> None:
    """Write code to files."""
    for source, suffix in zip(code, suffixes):
        with open(Path(output_dir) / (prefix + suffix), "w") as file:
            file.write(source)
