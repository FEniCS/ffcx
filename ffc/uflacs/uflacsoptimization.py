# -*- coding: utf-8 -*-
# Copyright (C) 2013-2017 Martin Sandve Alnæs
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import logging

logger = logging.getLogger(__name__)


def optimize_integral_ir(ir, parameters):
    """Compute optimized intermediate representation of integral."""

    logger.info("Optimizing uflacs representation")

    # TODO: Implement optimization of ssa representation prior to code
    # generation here.
    oir = ir

    return oir
