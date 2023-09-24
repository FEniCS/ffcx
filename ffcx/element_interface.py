# Copyright (C) 2021 Matthew W. Scroggs and Chris Richardson
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Finite element interface."""

import typing
import warnings

import basix
import basix.ufl
import ufl
import numpy as np
import numpy.typing as npt


def convert_element(element: ufl.finiteelement.FiniteElementBase) -> basix.ufl._ElementBase:
    """Convert and element to a FFCx element."""
    if isinstance(element, basix.ufl._ElementBase):
        return element
    else:
        warnings.warn(
            "Use of elements created by UFL is deprecated. You should create elements directly using Basix.",
            DeprecationWarning)
        return basix.ufl.convert_ufl_element(element)


def basix_index(indices: typing.Tuple[int]) -> int:
    """Get the Basix index of a derivative."""
    return basix.index(*indices)


def create_quadrature(
    cellname: str, degree: int, rule: str, elements: typing.List[basix.ufl._ElementBase]
) -> typing.Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Create a quadrature rule."""
    if cellname == "vertex":
        return (np.ones((1, 0), dtype=np.float64), np.ones(1, dtype=np.float64))
    else:
        celltype = basix.cell.string_to_type(cellname)
        polyset_type = basix.PolysetType.standard
        for e in elements:
            polyset_type = basix.polyset_superset(celltype, polyset_type, e.polyset_type)
        return basix.make_quadrature(
            celltype, degree, rule=basix.quadrature.string_to_type(rule), polyset_type=polyset_type)


def reference_cell_vertices(cellname: str) -> npt.NDArray[np.float64]:
    """Get the vertices of a reference cell."""
    return basix.geometry(basix.cell.string_to_type(cellname))


def map_facet_points(points: npt.NDArray[np.float64], facet: int, cellname: str) -> npt.NDArray[np.float64]:
    """Map points from a reference facet to a physical facet."""
    geom = basix.geometry(basix.cell.string_to_type(cellname))
    facet_vertices = [geom[i] for i in basix.topology(basix.cell.string_to_type(cellname))[-2][facet]]
    return np.asarray([facet_vertices[0] + sum((i - facet_vertices[0]) * j for i, j in zip(facet_vertices[1:], p))
                       for p in points], dtype=np.float64)


# TODO: remove this deprecated function
def QuadratureElement(
    cellname: str, value_shape: typing.Tuple[int, ...], scheme: typing.Optional[str] = None,
    degree: typing.Optional[int] = None, points: typing.Optional[npt.NDArray[np.float64]] = None,
    weights: typing.Optional[npt.NDArray[np.float64]] = None, mapname: str = "identity"
) -> basix.ufl._ElementBase:
    warnings.warn(
        "ffcx.element_interface.QuadratureElement is deprecated and will be removed after December 2023. "
        "Use basix.ufl.quadrature_element instead.", DeprecationWarning)
    return basix.ufl.quadrature_element(
        cell=cellname, value_shape=value_shape, scheme=scheme, degree=degree, points=points, weights=weights,
        mapname=mapname)


# TODO: remove this deprecated function
def RealElement(element: ufl.finiteelement.FiniteElementBase) -> basix.ufl._ElementBase:
    warnings.warn(
        "ffcx.element_interface.RealElement is deprecated and will be removed after December 2023. "
        "Use basix.ufl.real_element instead.", DeprecationWarning)
    return basix.ufl.real_element(cell=element.cell().cellname(), value_shape=element.value_shape())
