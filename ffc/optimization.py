# -*- coding: utf-8 -*-
# Copyright (C) 2009-2017 Anders Logg
#
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""
Compiler stage 5: optimization
------------------------------

This module implements the optimization of an intermediate code
representation.
"""
import logging

from ffc.representation import pick_representation

logger = logging.getLogger(__name__)


def optimize_ir(ir, parameters):
    "Optimize intermediate form representation."

    logger.info("Compiler stage 3: Optimizing intermediate representation")

    # Extract representations
    ir_elements, ir_dofmaps, ir_coordinate_mappings, ir_integrals, ir_forms = ir

    # Check if optimization is requested
    if not any(ir["integrals_metadata"]["optimize"] for ir in ir_integrals):
        logger.info(
            "Skipping optimizations, add -O or attach {{'optimize': True}} metadata to integrals")

    # Call on every bunch of integrals wich are compiled together
    oir_integrals = [
        _optimize_integral_ir(ir, parameters) if ir["integrals_metadata"]["optimize"] else ir
        for ir in ir_integrals
    ]

    return ir_elements, ir_dofmaps, ir_coordinate_mappings, oir_integrals, ir_forms


def _optimize_integral_ir(ir, parameters):
    "Compute optimized intermediate represention of integral."
    r = pick_representation(ir["representation"])
    return r.optimize_integral_ir(ir, parameters)
