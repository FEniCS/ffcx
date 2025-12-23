# Copyright (C) 2009-2025 Anders Logg, Martin Sandve Alnæs, Garth N. Wells, Chris Richardson and
# Paul T. Kühner
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# Note: Most of the code in this file is a direct translation from the
# old implementation in FFC
"""Generate file output for numba."""

import logging
import pprint
import textwrap

from numpy import typing as npt

from ffcx import __version__ as FFCX_VERSION
from ffcx.codegeneration import __version__ as UFC_VERSION
from ffcx.codegeneration.numba import file_template

logger = logging.getLogger("ffcx")


def generator(
    options: dict[str, int | float | npt.DTypeLike],
) -> tuple[tuple[str, str], tuple[str, str]]:
    """Generate UFCx code for file output.

    Args:
        options: Dict of options specified the kenerl generation, these will be documented in the
        generated file.

    Returns: tuple of file start- and end sections, each for declaration and implementation.

    """
    logger.info("Generating code for file")

    # Attributes
    d = {"ffcx_version": FFCX_VERSION, "ufcx_version": UFC_VERSION}
    d["options"] = textwrap.indent(pprint.pformat(options), "#  ")

    # Format declaration code
    code_pre = (
        file_template.declaration_pre.format_map(d),
        file_template.implementation_pre.format_map(d),
    )

    # Format implementation code
    code_post = (
        file_template.declaration_post.format_map(d),
        file_template.implementation_post.format_map(d),
    )

    return code_pre, code_post
