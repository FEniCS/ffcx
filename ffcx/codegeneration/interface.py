# Copyright (C) 2025 Paul T. Kühner
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

"""Backend interface declarations.

Every language backend needs to implement/overload this functionality.
"""

from collections.abc import Callable
from typing import Protocol

import basix
from numpy import typing as npt

import ffcx.codegeneration.lnodes as L
from ffcx.ir.representation import ExpressionIR, FormIR, IntegralIR


class Formatter(Protocol):
    """Formatter interface."""

    def __init__(self, dtype: npt.DTypeLike) -> None:
        """Create."""
        ...

    def __call__(self, obj: L.LNode) -> str:
        """Convert L-Node(s) to string representation."""
        ...


"""File to source string.

Note:
    Needs to be callable as file.generator.
"""
file_generator = Callable[[dict[str, int | float | npt.DTypeLike]], tuple[tuple[str], ...]]

"""Form to source string.

Note:
    Needs to be callable as form.generator.
"""
form_generator = Callable[[FormIR, dict[str, int | float | npt.DTypeLike]], tuple[str]]

"""Integral to source string.

Note:
    Needs to be callable as integral.generator.
"""
integral_generator = Callable[
    [IntegralIR, basix.CellType, dict[str, int | float | npt.DTypeLike]], tuple[str, ...]
]


"""Expression to source string.

Note:
    Needs to be callable as expression.generator.
"""
expression_generator = Callable[
    [ExpressionIR, dict[str, int | float | npt.DTypeLike]], tuple[str, ...]
]
