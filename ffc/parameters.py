# -*- coding: utf-8 -*-
# Copyright (C) 2005-2016 Anders Logg
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.

import os
import copy
from ffc.log import INFO


# Comments from other places in code:
# FIXME: Document option -fconvert_exceptions_to_warnings
# FIXME: Remove option epsilon and just rely on precision?

# NB! Parameters in the generate and build sets are
# included in jit signature, cache and log are not.
_FFC_GENERATE_PARAMETERS = {
    "format": "ufc",           # code generation format
    "representation": "auto",  # form representation / code
                               # generation strategy
    "quadrature_rule": "auto", # quadrature rule used for
                               # integration of element tensors
    "quadrature_degree": -1,   # quadrature degree used for
                               # computing integrals
    "precision": 0,            # precision used when writing
                               # numbers (0 for max precision)
    "epsilon": 1e-14,          # machine precision, used for
                               # dropping zero terms
                               # (quadrature and tensor repr)
    "split": False,            # split generated code into .h and
                               # .cpp file
    "form_postfix": True,      # postfix form name with "Function",
                               # "LinearForm" or BilinearForm
    "convert_exceptions_to_warnings": False,   # convert all exceptions to warning
                                               # in generated code
    "error_control": False,   # with error control
    "optimize": True,         # turn on optimization for code generation
    "max_signature_length": 0,  # set to positive integer to shorten signatures
    "generate_dummy_tabulate_tensor": False,  # set to True to replace tabulate_tensor body with no-op
    "add_tabulate_tensor_timing": False,      # set to True to add timing inside tabulate_tensor
}
_FFC_BUILD_PARAMETERS = {
    "cpp_optimize": True,          # optimization for the C++ compiler
    "cpp_optimize_flags": "-O2",   # optimization flags for the C++ compiler
}
_FFC_CACHE_PARAMETERS = {
    "cache_dir": "",        # cache dir used by Instant
    "output_dir": ".",      # output directory for generated code
}
_FFC_LOG_PARAMETERS = {
    "log_level": INFO + 5,  # log level, displaying only
                            # messages with level >= log_level
    "log_prefix": "",       # log prefix
}
FFC_PARAMETERS = {}
FFC_PARAMETERS.update(_FFC_BUILD_PARAMETERS)
FFC_PARAMETERS.update(_FFC_CACHE_PARAMETERS)
FFC_PARAMETERS.update(_FFC_LOG_PARAMETERS)
FFC_PARAMETERS.update(_FFC_GENERATE_PARAMETERS)


def split_parameters(parameters):
    """Split a parameters dict into groups based on what parameters are used for.

    """
    params = {
        "cache": {k: parameters[k] for k in _FFC_CACHE_PARAMETERS.keys()},
        "build": {k: parameters[k] for k in _FFC_BUILD_PARAMETERS.keys()},
        "generate": {k: parameters[k] for k in _FFC_GENERATE_PARAMETERS.keys()},
        "log": {k: parameters[k] for k in _FFC_LOG_PARAMETERS.keys()},
    }
    return params


def default_parameters():
    "Return (a copy of) the default parameter values for FFC."
    return copy.deepcopy(FFC_PARAMETERS)


def default_jit_parameters():
    parameters = default_parameters()

    # TODO: This is not in the above parameters dict.
    #       There are other parameters like this.
    #       This is confusing, which parameters are available? What are the defaults?
    # Skip evaluation of basis derivatives in elements by default because it's costly
    # FIXME: Make this False when we have elements generated once instead of for each form
    parameters["no-evaluate_basis_derivatives"] = True

    # Don't postfix form names
    parameters["form_postfix"] = False

    return parameters


def validate_parameters(parameters):
    "Initial check of parameters."
    p = default_parameters()
    if parameters is not None:
        p.update(parameters)

    #if isinstance(p["optimize"], bool):
    #    p["optimize"] = 1 if p["optimize"] else 0
    if isinstance(p["optimize"], int):
        p["optimize"] = bool(p["optimize"])

    return p


def validate_jit_parameters(parameters):
    "Check parameters and add any missing parameters"
    p = default_jit_parameters()
    if parameters is not None:
        p.update(parameters)

    #if isinstance(p["optimize"], bool):
    #    p["optimize"] = 1 if p["optimize"] else 0
    if isinstance(p["optimize"], int):
        p["optimize"] = bool(p["optimize"])

    return p


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


def compute_jit_parameters_signature(parameters):
    "Return parameters signature (some parameters must be ignored)."
    from ufl.utils.sorting import canonicalize_metadata
    parameters = compilation_relevant_parameters(parameters)
    return str(canonicalize_metadata(parameters))

