# Copyright (C) 2016 Jan Blechta
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import logging

from ffc.ir.representationutils import initialize_integral_ir

logger = logging.getLogger(__name__)


def compute_integral_ir(integral_data, form_data, form_id, element_numbers, classnames, parameters):
    """Compute intermediate represention of integral."""

    logger.info("Computing tsfc representation")

    # Initialise representation
    ir = initialize_integral_ir("tsfc", integral_data, form_data, form_id)

    # TSFC treats None and unset differently, so remove None values.
    parameters = {k: v for k, v in parameters.items() if v is not None}

    # Delay TSFC compilation
    ir["compile_integral"] = (integral_data, form_data, None, parameters)

    return ir
