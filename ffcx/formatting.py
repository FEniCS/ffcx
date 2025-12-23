# Copyright (C) 2009-2018 Anders Logg and Garth N. Wells
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Compiler stage 5: Code formatting.

This module implements the formatting of UFCx code from a given
dictionary of generated C++ code for the body of each UFCx function.

It relies on templates for UFCx code available as part of the module
ufcx_utils.
"""

import logging
from pathlib import Path

from ffcx.codegeneration.codegeneration import CodeBlocks

logger = logging.getLogger("ffcx")


def format_code(code: CodeBlocks) -> tuple[str, str]:
    """Format given code in UFCx format.

    Returns two strings with header and source file contents.
    """
    logger.info(79 * "*")
    logger.info("Compiler stage 5: Formatting code")
    logger.info(79 * "*")

    code_c = ""
    code_h = ""
    for parts_code in code:
        code_h += "".join([c[0] for c in parts_code])
        code_c += "".join([c[1] for c in parts_code])

    return code_h, code_c


def write_code(
    code_h: str,
    code_c: str,
    filename_stem: str,
    suffixes: tuple[str | None, str | None],
    output_dir: str,
) -> None:
    """Write code to files.

    Args:
        code_h: Header file content.
        code_c: Source file content.
        filename_stem: The stem of the filename to use for both header
            and source files.
        suffixes: Declaration/implementation file suffixes.
        output_dir: Directory where the files should be written.
    """
    # def _write_file(output: str, prefix: str, suffix: str, output_dir: str) -> None:
    #     """Write generated code to file."""
    #     with open(Path(output_dir) / (prefix + suffix), "w") as file:
    #         file.write(output)

    def _write_file(output: str, prefix: str, postfix: str, output_dir: str) -> None:
        """Write generated code to file."""
        filename = Path(output_dir, prefix).with_suffix(postfix)
        assert filename.parent.exists(), f"Output directory '{filename.parent}' does not exist."
        filename.write_text(output)

    if suffixes[0] is not None:
        _write_file(code_h, filename_stem, suffixes[0], output_dir)
    if suffixes[1] is not None:
        _write_file(code_c, filename_stem, suffixes[1], output_dir)
