# Copyright (C) 2025 Paul T. Kühner
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

"""Backend interface declarations."""

import abc

from numpy import typing as npt

import ffcx.codegeneration.lnodes as L


class Formatter(abc.ABC):
    """Formatter interface."""

    def __init__(self, dtype: npt.DTypeLike) -> None:
        """Create."""
        raise NotImplementedError

    def __call__(self, obj: L.LNode) -> str:
        """Convert L-Node(s) to string representation."""
        raise NotImplementedError
