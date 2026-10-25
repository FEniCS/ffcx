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
    "language": (str, "C", "Language to output kernel in", None),
    "scalar_type": (
        str,
        "float64",
        "Scalar type to use in generated code.",
        ("float32", "float64", "complex64", "complex128"),
    ),
    "sum_factorization": (
        bool,
        False,
        "Use sum factorization for cell integrals of tensor product elements: tabulate on a "
        "tensor product quadrature rule and contract one direction at a time, at a cost of "
        "O(P^(2d+1)) for a bilinear form of degree P in dimension d.",
        None,
    ),
    "table_rtol": (
        float,
        1e-6,
        "Relative precision to use when comparing finite element table values for reuse.",
        None,
    ),
    "table_atol": (
        float,
        1e-9,
        "Absolute precision to use when comparing finite element table values reuse.",
        None,
    ),
    "verbosity": (
        int,
        30,
        "Logger verbosity, follows standard library levels, i.e. INFO=20, DEBUG=10, etc.",
        None,
    ),
    "part": (str, "full", "Part of bilinear tensor to assemble", ("full", "diagonal")),
    "do_cancel_jacobian_products": (
        bool,
        True,
        "Cancel contractions of the Jacobian with its inverse in form integrands.",
        None,
    ),
}


def _check_option_keys(options: dict, source: str) -> None:
    """Check that all keys are recognised FFCx options.

    Args:
        options: Option values to check.
        source: Description of where the options came from, used in the
            error message.

    Raises:
        ValueError: If any key is not in `FFCX_DEFAULT_OPTIONS`.
    """
    unknown = sorted(set(options) - set(FFCX_DEFAULT_OPTIONS))
    if unknown:
        known = ", ".join(sorted(FFCX_DEFAULT_OPTIONS))
        raise ValueError(
            f"Unknown FFCx option(s) in {source}: {', '.join(unknown)}. "
            f"Available options are: {known}."
        )


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

        A key that is not in `FFCX_DEFAULT_OPTIONS` raises `ValueError`,
        from whichever of these sources it came.

        Example `ffcx_options.json` file:

          { "table_atol": 1e-10 }

    """
    options: dict[str, npt.DTypeLike | int | float] = {}

    for opt, (_, value, _, _) in FFCX_DEFAULT_OPTIONS.items():
        options[opt] = value  # type: ignore

    # NOTE: _load_options uses functools.lru_cache
    user_options, pwd_options = _load_options()

    # Reject unrecognised keys. An unknown key would otherwise be merged
    # in silently: it has no effect on the generated code, but it does
    # change the JIT cache signature, so the form is recompiled and the
    # option appears to have been applied.
    _check_option_keys(user_options, "$XDG_CONFIG_HOME/ffcx/ffcx_options.json")
    _check_option_keys(pwd_options, "$PWD/ffcx_options.json")
    if priority_options is not None:
        _check_option_keys(priority_options, "the options passed to FFCx")

    options.update(user_options)
    options.update(pwd_options)
    if priority_options is not None:
        options.update(priority_options)

    logger.setLevel(int(options["verbosity"]))  # type: ignore
    # pformat dominates this function and runs on every JIT cache hit.
    if logger.isEnabledFor(logging.INFO):
        logger.info("Final option values")
        logger.info(pprint.pformat(options))

    return options


def get_language(options: dict[str, int | float | npt.DTypeLike]) -> str:
    """Retrieve the language option from the options database.

    Applies for internal languages the alias conversion.
    """
    lang = str(options.get("language", FFCX_DEFAULT_OPTIONS["language"][1]))
    return {"C": "ffcx.codegeneration.C", "numba": "ffcx.codegeneration.numba"}.get(lang, lang)
