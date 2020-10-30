# Copyright (C) 2005-2020 Anders Logg, Michal Habera, Jack S. Hale
#
# This file is part of FFCX. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import json
import logging
import os
import os.path
import pathlib

logger = logging.getLogger("ffcx")

FFCX_DEFAULT_PARAMETERS = {
    "quadrature_rule":
        ("default", "Quadrature rule used for integration of element tensors."),
    "quadrature_degree":
        (-1, """Quadrature degree used for approximating integrals.
                (-1 means automatically determined from integrand polynomial degree)"""),
    "precision":
        (-1, """Precision used when writing numbers (-1 for max precision).
                Represents maximum number of digits after decimal point."""),
    "epsilon":
        (1e-14, "Machine precision, used for dropping zero terms in tables"),
    "scalar_type":
        ("double", "Scalar type used in generated code. Any of real or complex C floating-point types."),
    "tabulate_tensor_void":
        (False, "True to generate empty tabulation kernels."),
    "table_rtol":
        (1e-6, "Relative precision to use when comparing finite element table values for table reuse."),
    "table_atol":
        (1e-9, "Absolute precision to use when comparing finite element table values for reuse."),
    "assume_aligned":
        (-1, """Assumes alignment (in bytes) of pointers to tabulated tensor, coefficients and constants array.
               This value must be compatible with alignment of data structures allocated outside FFC.
               (-1 means no alignment assumed, safe option)"""),
    "padlen":
        (1, "Pads every declared array in tabulation kernel such that its last dimension is divisible by given value."),
    "verbosity":
        (30, "Logger verbosity. Follows standard logging library levels, i.e. INFO=20, DEBUG=10, etc.")
}


def get_parameters():
    """Return (a copy of) the parameter values for FFCX.

    Priority order is:
      environment variables (to discuss)
      ~/.config/ffcx/paramters.json (user parameters)
      $(pwd)/.ffcx_parameters.json (local parameters)
      FFCX_DEFAULT_PARAMETERS
    """
    parameters = {}

    for param, (value, desc) in FFCX_DEFAULT_PARAMETERS.items():
        parameters[param] = value

    user_config_file = os.path.join(pathlib.Path.home(), ".config", "ffcx", "parameters.json")
    pwd_config_file = os.path.join(os.getcwd(), ".ffcx_parameters.json")

    try:
        with open(user_config_file) as f:
            user_parameters = json.load(f)
    except FileNotFoundError:
        user_parameters = {}

    try:
        with open(pwd_config_file) as f:
            pwd_parameters = json.load(f)
    except FileNotFoundError:
        pwd_parameters = {}

    keys = os.environ.keys()
    env_parameters = {}
    for name, value in FFCX_DEFAULT_PARAMETERS.items():
        param_name = "FFCX_" + name.upper()
        param_type = type(value[0])

        if param_name in keys:
            env_parameters[name] = param_type(os.environ[param_name])
            logger.info("Parameter {} forced to {} from environmental variable.".format(name, env_parameters[name]))

    parameters.update(user_parameters)
    parameters.update(pwd_parameters)
    parameters.update(env_parameters)

    logging.debug("Final parameter settings")
    logging.debug(parameters)

    return parameters
