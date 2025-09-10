# Copyright (C) 2025 Jørgen S. Dokken
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

"""Module for storing type definitions used in the FFCx code base."""

from typing import NamedTuple


class IntegralType(NamedTuple):
    """Class for storing information about an integral type."""

    codim: int  # Codimension of integral with respect to the topological dimension of the domain
    num_cells: int = 1  # Number of cells the integral is over
    is_expression: bool = (
        False  # True if the integral is an Expression (integrand) rather than an integral
    )
    is_custom: bool = False  # True if the integral uses the custom ufl integral measure.


def convert_to_integral_type(integral_type: str) -> IntegralType:
    """Convert UFL integral type string to IntegralType."""
    match integral_type:
        case "cell":
            return IntegralType(codim=0, num_cells=1)
        case "exterior_facet":
            return IntegralType(codim=1, num_cells=1)
        case "interior_facet":
            return IntegralType(codim=1, num_cells=2)
        case "ridge":
            return IntegralType(codim=2, num_cells=1)
        case "vertex":
            return IntegralType(codim=-1, num_cells=1)
        case "custom":
            return IntegralType(codim=0, num_cells=1, is_custom=True)
        case _:
            raise ValueError(f"Unknown integral type:{integral_type}")
