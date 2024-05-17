# Copyright (C) 2005-2020 Anders Logg, Michal Habera, Jack S. Hale
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Options."""

from __future__ import annotations

import functools
import json
import logging
import os
import os.path
import pprint
from pathlib import Path

import numpy.typing as npt

logger = logging.getLogger("ffcx")

FFCX_DEFAULT_OPTIONS = {
    "epsilon": (float, 1e-14, "machine precision, used for dropping zero terms in tables.", None),
    "scalar_type": (
        str,
        "float64",
        "scalar type to use in generated code.",
        ("float32", "float64", "complex64", "complex128"),
    ),
    "sum_factorization": (bool, False, "use sum factorization.", None),
    "table_rtol": (
        float,
        1e-6,
        "relative precision to use when comparing finite element table values for reuse.",
        None,
    ),
    "table_atol": (
        float,
        1e-9,
        "absolute precision to use when comparing finite element table values reuse.",
        None,
    ),
    "verbosity": (
        int,
        30,
        "logger verbosity, follows standard library levels, i.e. INFO=20, DEBUG=10, etc.",
        None,
    ),
}


@functools.cache
def _load_options() -> tuple[dict, dict]:
    """Load options from JSON files."""
    user_config_file = os.getenv("XDG_CONFIG_HOME", default=Path.home().joinpath(".config")) / Path(
        "ffcx", "ffcx_options.json"
    )
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


def get_options(
    priority_options: dict[str, npt.DTypeLike | int | float] | None = None,
) -> dict[str, int | float | npt.DTypeLike]:
    """Return (a copy of) the merged option values for FFCX.

    Args:
        priority_options: take priority over all other option values (see notes)

    Returns:
        merged option values

    Note:
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

          { "epsilon": 1e-7 }

    """
    options: dict[str, npt.DTypeLike | int | float] = {}

    for opt, (_, value, _, _) in FFCX_DEFAULT_OPTIONS.items():
        options[opt] = value  # type: ignore

    # NOTE: _load_options uses functools.lru_cache
    user_options, pwd_options = _load_options()

    options.update(user_options)
    options.update(pwd_options)
    if priority_options is not None:
        options.update(priority_options)

    logger.setLevel(int(options["verbosity"]))  # type: ignore
    logger.info("Final option values")
    logger.info(pprint.pformat(options))

    return options
