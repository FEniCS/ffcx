# Copyright (C) 2005-2020 Anders Logg, Michal Habera
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import copy
import logging

logger = logging.getLogger("ffcx")

FFCX_PARAMETERS = {
    "quadrature_rule": "default",  # quadrature rule used for integration of element tensors
    "quadrature_degree": -1,  # quadrature degree used for computing integrals (-1 means auto)
    "precision": -1,  # precision used when writing numbers (-1 for max precision)
    "epsilon": 1e-14,  # machine precision, used for dropping zero terms in tables
    # Scalar type to be used in generated code (real or complex
    # C double precision floating-point types)
    "scalar_type": "double",
    "external_includes": "",  # ':' separated list of include filenames to add to generated code
    "tabulate_tensor_void": False,  # generate empty tabulation kernels, for benchmarking

    # Relative precision to use when comparing finite element table
    # values for table reuse
    "table_rtol": 1e-6,

    # Absolute precision to use when comparing finite element table
    # values for table reuse and dropping of table zeros
    "table_atol": 1e-9,

    # Number of points to evaluate
    "chunk_size": 8,

    "alignas": 32,
    "assume_aligned": None,
    "padlen": 1,
    "use_symbol_array": True
}


def default_parameters():
    """Return (a copy of) the default parameter values for FFCX."""
    parameters = copy.deepcopy(FFCX_PARAMETERS)

    return parameters
