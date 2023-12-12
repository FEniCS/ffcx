# Copyright (C) 2009-2022 Anders Logg, Martin Sandve Alnæs, Matthew Scroggs
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

# Note: Much of the code in this file is a direct translation
# from the old implementation in FFC, although some improvements
# have been made to the generated code.

import logging

# import ffcx.codegeneration.C.basix_custom_element_template as ufcx_basix_custom_finite_element
# import ffcx.codegeneration.C.finite_element_template as ufcx_finite_element
# import ufl

logger = logging.getLogger("ffcx")


def generator(ir, options):
    """Generate UFC code for a finite element."""
    logger.info("Generating code for finite element:")
    logger.info(f"--- degree: {ir.degree}")
    logger.info(f"--- value shape: {ir.value_shape}")
    logger.info(f"--- name: {ir.name}")

    return ("", "")
