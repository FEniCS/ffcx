# Copyright (C) 2025 Paul T. Kühner
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

"""Backend interface declarations."""

from typing import Protocol

import basix
from numpy import typing as npt

import ffcx.codegeneration.lnodes as L
from ffcx.ir.representation import FormIR, IntegralIR


class Formatter(Protocol):
    """Formatter interface."""

    def __init__(self, dtype: npt.DTypeLike) -> None:
        """Create."""
        ...

    def __call__(self, obj: L.LNode) -> str:
        """Convert L-Node(s) to string representation."""
        ...


class file(Protocol):
    """File formatter."""

    def generator(options: dict[str, int | float | npt.DTypeLike]) -> tuple[str, ...]:
        """File to source string."""
        ...


class form(Protocol):
    """Form formatter."""

    def generator(ir: FormIR, options: dict[str, int | float | npt.DTypeLike]) -> tuple[str]:
        """Form to source string."""
        ...


class integral(Protocol):
    """Integral formatter."""

    def generator(
        ir: IntegralIR, domain: basix.CellType, options: dict[str, int | float | npt.DTypeLike]
    ) -> tuple[str]:
        """Integral to source string."""
        ...
