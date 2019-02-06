# -*- coding: utf-8 -*-
# Copyright (C) 2005-2017 Anders Logg
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import copy
import logging
import os

logger = logging.getLogger(__name__)

# Comments from other places in code:
# FIXME: Document option -fconvert_exceptions_to_warnings

# NB! Parameters in the generate and build sets are
# included in jit signature, cache and log are not.
_FFC_GENERATE_PARAMETERS = {
    "format": "ufc",  # code generation format
    "representation": "auto",  # form representation / code generation strategy
    "quadrature_rule":
    # quadrature rule used for integration of element tensors (None is auto)
    None,
    # quadrature degree used for computing integrals (None is auto)
    "quadrature_degree": None,
    # precision used when writing numbers (None for max precision)
    "precision": None,
    "epsilon": 1e-14,  # machine precision, used for dropping zero terms in tables
    "form_postfix": True,  # postfix form name with "Function", "LinearForm" or BilinearForm
    # convert all exceptions to warning in generated code
    "convert_exceptions_to_warnings": False,
    "max_signature_length":
    0,  # set to positive integer to shorten signatures set to True to replace tabulate_tensor body with no-op
    # set to True to add timing inside tabulate_tensor
    "generate_dummy_tabulate_tensor": False,
    "add_tabulate_tensor_timing": False,
    # Scalar type to be used in generated code (real or complex
    # C double precision floating-point types)
    "scalar_type": "double",
    # Max time to wait on cache if not building on this
    # process (seconds)
    "timeout": 10,
    # ':' separated list of include filenames to add to generated code
    "external_includes": "",
    # Whether to crosslink JIT libraries or build standalone
    "crosslink": False,
}
_FFC_BUILD_PARAMETERS = {
    "cpp_optimize": True,  # optimization for the C++ compiler
    "cpp_optimize_flags": "-O2",  # optimization flags for the C++ compiler
    # ':' separated list of libraries to link JIT compiled libraries with
    "external_libraries": "",
    "external_library_dirs":
    "",  # ':' separated list of library search dirs to add when JIT compiling
    # ':' separated list of include dirs to add when JIT compiling
    "external_include_dirs": "",
}
_FFC_CACHE_PARAMETERS = {
    "cache_dir": "compile_cache",  # cache dir used by default
    "output_dir": ".",  # output directory for generated code
}
_FFC_LOG_PARAMETERS = {
    # "log_level": INFO + 5,  # log level, displaying only messages with level >= log_level
    "log_prefix": "",  # log prefix
    "visualise": False,
}
FFC_PARAMETERS = {}
FFC_PARAMETERS.update(_FFC_BUILD_PARAMETERS)
FFC_PARAMETERS.update(_FFC_CACHE_PARAMETERS)
FFC_PARAMETERS.update(_FFC_LOG_PARAMETERS)
FFC_PARAMETERS.update(_FFC_GENERATE_PARAMETERS)


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

    # Don't postfix form names
    parameters["form_postfix"] = False

    return parameters


def validate_parameters(parameters):
    """Initial check of parameters."""
    p = default_parameters()
    if parameters is not None:
        p.update(parameters)

    _validate_parameters(p)

    return p


def validate_jit_parameters(parameters):
    """Check parameters and add any missing parameters"""
    p = default_jit_parameters()
    if parameters is not None:
        p.update(parameters)

    _validate_parameters(p)

    return p


def _validate_parameters(parameters):
    """Does some casting of parameter values in place on the
    provided dictionary"""
    # Convert all legal default values to None
    if parameters["quadrature_rule"] in ["auto", None, "None"]:
        parameters["quadrature_rule"] = None

    # Convert all legal default values to None and
    # cast nondefaults from str to int
    if parameters["quadrature_degree"] in ["auto", -1, None, "None"]:
        parameters["quadrature_degree"] = None
    else:
        try:
            parameters["quadrature_degree"] = int(parameters["quadrature_degree"])
        except Exception:
            logger.exception("Failed to convert quadrature degree '%s' to int" %
                             parameters.get("quadrature_degree"))
            raise

    # Convert all legal default values to None and
    # cast nondefaults from str to int
    if parameters["precision"] in ["auto", None, "None"]:
        parameters["precision"] = None
    else:
        try:
            parameters["precision"] = int(parameters["precision"])
        except Exception:
            logger.exception("Failed to convert precision '{}' to int".format(
                parameters.get("precision")))
            raise


def compilation_relevant_parameters(parameters):
    p = parameters.copy()
    for k in _FFC_LOG_PARAMETERS:
        del p[k]
    for k in _FFC_CACHE_PARAMETERS:
        del p[k]

    # This doesn't work because some parameters may not be among the defaults above.
    # That is somewhat confusing but we'll just have to live with it at least for now.
    # sp = split_parameters(parameters)
    # p = {}
    # p.update(sp["generate"])
    # p.update(sp["build"])

    return p


def compute_jit_signature(parameters):
    """Return parameters signature (some parameters must be ignored)."""
    from ufl.utils.sorting import canonicalize_metadata
    parameters = compilation_relevant_parameters(parameters)
    return str(canonicalize_metadata(parameters))
