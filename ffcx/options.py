# Copyright (C) 2005-2020 Anders Logg, Michal Habera, Jack S. Hale
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import functools
import json
import logging
import os
import os.path
import pprint
from typing import Optional, Dict, Any
from pathlib import Path

logger = logging.getLogger("ffcx")

FFCX_DEFAULT_OPTIONS = {
    "epsilon":
        (1e-14, "Machine precision, used for dropping zero terms in tables"),
    "scalar_type":
        ("double", """Scalar type used in generated code. Any of real or complex C floating-point types, e.g.
                      float, double, float _Complex, double _Complex, ..."""),
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


@functools.lru_cache(maxsize=None)
def _load_options():
    """Load options from JSON files."""
    user_config_file = os.getenv("XDG_CONFIG_HOME", default=Path.home().joinpath(".config")) \
        / Path("ffcx", "ffcx_options.json")
    try:
        with open(user_config_file) as f:
            user_options = json.load(f)
    except FileNotFoundError:
        user_options = {}

    pwd_config_file = Path.cwd().joinpath("ffcx_options.json")
    try:
        with open(pwd_config_file) as f:
            pwd_options = json.load(f)
    except FileNotFoundError:
        pwd_options = {}

    return (user_options, pwd_options)


def get_options(priority_options: Optional[dict] = None) -> dict:
    """Return (a copy of) the merged option values for FFCX.

    Options
    ----------
      priority_options:
        take priority over all other option values (see notes)

    Returns
    -------
      dict: merged option values

    Notes
    -----
    This function sets the log level from the merged option values prior to
    returning.

    The `ffcx_options.json` files are cached on the first call. Subsequent
    calls to this function use this cache.

    Priority ordering of options from highest to lowest is:

    -  **priority_options** (API and command line options)
    -  **$PWD/ffcx_options.json** (local options)
    -  **$XDG_CONFIG_HOME/ffcx/ffcx_options.json** (user options)
    -  **FFCX_DEFAULT_OPTIONS** in `ffcx.options`

    `XDG_CONFIG_HOME` is `~/.config/` if the environment variable is not set.

    Example `ffcx_options.json` file:

      { "assume_aligned": 32, "epsilon": 1e-7 }

    """
    options: Dict[str, Any] = {}

    for opt, (value, _) in FFCX_DEFAULT_OPTIONS.items():
        options[opt] = value

    # NOTE: _load_options uses functools.lru_cache
    user_options, pwd_options = _load_options()

    options.update(user_options)
    options.update(pwd_options)
    if priority_options is not None:
        options.update(priority_options)

    logger.setLevel(options["verbosity"])

    logger.info("Final option values")
    logger.info(pprint.pformat(options))

    return options
