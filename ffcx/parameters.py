# Copyright (C) 2005-2020 Anders Logg, Michal Habera
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import copy
import logging
import os

logger = logging.getLogger(__name__)

FFCX_PARAMETERS = {
    "representation": None,  # form representation / code generation strategy
    "quadrature_rule": "auto",  # quadrature rule used for integration of element tensors (None is auto)
    "quadrature_degree": "auto",  # quadrature degree used for computing integrals (None is auto)
    "precision": "max",  # precision used when writing numbers (None for max precision)
    "epsilon": 1e-14,  # machine precision, used for dropping zero terms in tables
    # Scalar type to be used in generated code (real or complex
    # C double precision floating-point types)
    "scalar_type": "double",
    "external_includes": "",  # ':' separated list of include filenames to add to generated code
}


def default_parameters():
    """Return (a copy of) the default parameter values for FFCX."""
    parameters = copy.deepcopy(FFCX_PARAMETERS)

    return parameters
