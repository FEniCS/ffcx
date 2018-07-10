# -*- coding: utf-8 -*-
# Copyright (C) 2013-2017 Martin Sandve Alnæs
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Controlling algorithm for building the tabulate_tensor
source structure from factorized representation."""

import logging

from ffc.backends.ffc.backend import FFCBackend
from ffc.representationutils import initialize_integral_code
from ffc.uflacs.integralgenerator import IntegralGenerator
from ffc.uflacs.language.format_lines import format_indented_lines

logger = logging.getLogger(__name__)


def generate_integral_code(ir, prefix, parameters):
    """Generate code for integral from intermediate representation."""

    logger.info("Generating code from ffc.uflacs representation")

    # FIXME: Is this the right precision value to use? Make it default to None or 0.
    precision = ir["integrals_metadata"]["precision"]

    # Create FFC C backend
    backend = FFCBackend(ir, parameters)

    # Configure kernel generator
    ig = IntegralGenerator(ir, backend, precision)

    # Generate code ast for the tabulate_tensor body
    parts = ig.generate()

    # Format code as string
    body = format_indented_lines(parts.cs_format(precision), 1)

    # Generate generic ffc code snippets and add uflacs specific parts
    code = initialize_integral_code(ir, prefix, parameters)
    code["tabulate_tensor"] = body
    code["additional_includes_set"] = set(ig.get_includes())
    code["additional_includes_set"].update(ir.get("additional_includes_set", ()))

    return code
