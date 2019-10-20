# Copyright (C) 2005-2017 Anders Logg
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import copy
import logging
import os

logger = logging.getLogger(__name__)

FFC_PARAMETERS = {
    "representation": "auto",  # form representation / code generation strategy
    "quadrature_rule": None,  # quadrature rule used for integration of element tensors (None is auto)
    "quadrature_degree": None,  # quadrature degree used for computing integrals (None is auto)
    "precision": None,  # precision used when writing numbers (None for max precision)
    "epsilon": 1e-14,  # machine precision, used for dropping zero terms in tables
    # Scalar type to be used in generated code (real or complex
    # C double precision floating-point types)
    "scalar_type": "double",
    "external_includes": "",  # ':' separated list of include filenames to add to generated code
}


def default_parameters():
    """Return (a copy of) the default parameter values for FFC."""
    parameters = copy.deepcopy(FFC_PARAMETERS)

    # HACK
    r = os.environ.get("FFC_FORCE_REPRESENTATION")
    if r:
        parameters["representation"] = r

    return parameters


def default_jit_parameters():
    parameters = default_parameters()
    return parameters


def validate_parameters(parameters):
    """Initial check of parameters."""
    p = default_parameters()
    if parameters is not None:
        p.update(parameters)
    _validate_parameters(p)
    return p


def _validate_parameters(parameters):
    """Does some casting of parameter values in place on the provided dictionary."""
    # Convert all legal default values to None
    if parameters["quadrature_rule"] in ("auto", None, "None"):
        parameters["quadrature_rule"] = None

    # Convert all legal default values to None and cast nondefaults from
    # str to int
    if parameters["quadrature_degree"] in ("auto", -1, None, "None"):
        parameters["quadrature_degree"] = None
    else:
        try:
            parameters["quadrature_degree"] = int(parameters["quadrature_degree"])
        except Exception:
            logger.exception("Failed to convert quadrature degree '%s' to int" %
                             parameters.get("quadrature_degree"))
            raise

    # Convert all legal default values to None and cast nondefaults from
    # str to int
    if parameters["precision"] in ["auto", None, "None"]:
        parameters["precision"] = None
    else:
        try:
            parameters["precision"] = int(parameters["precision"])
        except Exception:
            logger.exception("Failed to convert precision '{}' to int".format(
                parameters.get("precision")))
            raise
