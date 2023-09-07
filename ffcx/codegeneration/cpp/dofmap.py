# Copyright (C) 2009-2018 Anders Logg, Martin Sandve Alnæs and Garth N. Wells
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# Note: Most of the code in this file is a direct translation from the
# old implementation in FFC

import logging

# import ffcx.codegeneration.cpp.dofmap_template as ufcx_dofmap

logger = logging.getLogger("ffcx")


def generator(ir, options):
    """Generate UFC code for a dofmap."""
    logger.info("Generating code for dofmap:")
    logger.info(f"--- num element support dofs: {ir.num_element_support_dofs}")
    logger.info(f"--- name: {ir.name}")

    return ("", "")
