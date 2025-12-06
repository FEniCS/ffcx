# Copyright (C) 2009-2017 Anders Logg and Martin Sandve Alnæs
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

# Note: Most of the code in this file is a direct translation from the
# old implementation in FFC
"""Generate code for a form."""

import logging

# from ffcx.codegeneration.cpp import form_template

logger = logging.getLogger("ffcx")


def generator(ir, options):
    """Generate UFC code for a form."""
    logger.info("Generating code for form:")
    logger.info(f"--- rank: {ir.rank}")
    logger.info(f"--- name: {ir.name}")

    return ("", "")
