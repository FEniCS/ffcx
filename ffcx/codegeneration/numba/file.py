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
from ffcx.codegeneration.common import template_keys
from ffcx.codegeneration.numba import file_template

suffixes = ("_numba.py",)

logger = logging.getLogger("ffcx")


def generator(
    options: dict[str, int | float | npt.DTypeLike],
) -> tuple[tuple[str], tuple[str]]:
    """Generate UFC code for file output.

    Args:
        options: Dict of options specified the kenerl generation, these will be documented in the
        generated file.

    Returns: tuple of file start- and end sections, each for declaration and implementation.

    """
    logger.info("Generating code for file")

    # Attributes
    d = {"ffcx_version": FFCX_VERSION, "ufcx_version": UFC_VERSION}
    d["options"] = textwrap.indent(pprint.pformat(options), "#  ")

    assert set(d.keys()) == template_keys(file_template.factory)
    return (file_template.factory.format_map(d),), ("",)
