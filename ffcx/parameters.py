# Copyright (C) 2005-2020 Anders Logg, Michal Habera, Jack S. Hale
#
# This file is part of FFCX. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import configparser
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
        (1, "Pads every declared array in tabulation kernel such that its last dimension is divisible by given value.")
}


def get_parameters():
    """Return (a copy of) the parameter values for FFCX.

    Priority order is:
      Environment variables
      ~/.config/ffcx/paramters.cfg (user parameters)
      $(pwd)/.ffcx_parameters.cfg (local parameters)
      FFCX_DEFAULT_PARAMETERS

    according to the logic of the standard Python ConfigParser class.
    """
    parameters = {}

    for param, (value, desc) in FFCX_DEFAULT_PARAMETERS.items():
        parameters[param] = value

    USER_CONFIG_FILE = os.path.join(pathlib.Path.home(), ".config", "ffcx", "parameters.cfg")
    PWD_CONFIG_FILE = os.path.join(os.getcwd(), ".ffcx_parameters.cfg")

    config = configparser.ConfigParser()
    parsed_config_files = config.read([parameters, USER_CONFIG_FILE, PWD_CONFIG_FILE])

    for parsed_config_file in parsed_config_files:
        logging.debug("Using parameter file at: ".format(str(parsed_config_file)))

    logging.debug("Final parameter settings")
    logging.debug(parameters)

    return parameters
