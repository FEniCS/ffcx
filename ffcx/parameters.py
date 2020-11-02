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
import pprint
from typing import Optional

logger = logging.getLogger("ffcx")

FFCX_DEFAULT_PARAMETERS = {
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


def get_parameters(priority_parameters: Optional[dict] = None) -> dict:
    """Return (a copy of) the merged parameter values for FFCX.

    Parameters
    ----------
      priority_parameters: 
        take priority over all other parameter values (see notes)

    Returns
    -------
      dict: merged parameter values

    Notes
    -----

    This function sets the log level from the merged parameter values prior to
    returning.

    Priority ordering of parameters is:
      priority_parameters (API and command line parameters)
      $(pwd)/ffcx_parameters.json (local parameters)
      ~/.config/ffcx/ffcx_paramters.json (user parameters)
      FFCX_DEFAULT_PARAMETERS in ffcx.parameters

    Example ffcx_parameters.json file:

      { "assume_aligned": 32, "epsilon": 1e-7 }

    """
    parameters = {}

    for param, (value, desc) in FFCX_DEFAULT_PARAMETERS.items():
        parameters[param] = value

    user_config_file = os.path.join(pathlib.Path.home(), ".config", "ffcx", "ffcx_parameters.json")
    pwd_config_file = os.path.join(os.getcwd(), "ffcx_parameters.json")

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

    parameters.update(user_parameters)
    parameters.update(pwd_parameters)
    parameters.update(priority_parameters)

    logger.setLevel(parameters["verbosity"])

    logger.info("Final parameter values")
    logger.info(pprint.pformat(parameters))

    return parameters
