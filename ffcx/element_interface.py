# Copyright (C) 2021 Matthew W. Scroggs and Chris Richardson
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Finite element interface."""

import basix
import basix.ufl
import numpy as np
import numpy.typing as npt
from basix import CellType as _CellType


def basix_index(indices: tuple[int]) -> int:
    """Get the Basix index of a derivative."""
    return basix.index(*indices)


def create_quadrature(
    cellname: str, degree: int, rule: str, elements: list[basix.ufl._ElementBase]
) -> tuple[npt.ArrayLike, npt.ArrayLike]:
    """Create a quadrature rule."""
    if cellname == "vertex":
        return (np.ones((1, 0), dtype=np.float64), np.ones(1, dtype=np.float64))
    else:
        celltype = _CellType[cellname]
        polyset_type = basix.PolysetType.standard
        for e in elements:
            polyset_type = basix.polyset_superset(celltype, polyset_type, e.polyset_type)
        return basix.make_quadrature(
            celltype, degree, rule=basix.quadrature.string_to_type(rule), polyset_type=polyset_type
        )


def reference_cell_vertices(cellname: str) -> npt.NDArray[np.float64]:
    """Get the vertices of a reference cell."""
    return np.asarray(basix.geometry(_CellType[cellname]))


def map_facet_points(
    points: npt.NDArray[np.float64], facet: int, cellname: str
) -> npt.NDArray[np.float64]:
    """Map points from a reference facet to a physical facet."""
    geom = np.asarray(basix.geometry(_CellType[cellname]))
    facet_vertices = [geom[i] for i in basix.topology(_CellType[cellname])[-2][facet]]
    return np.asarray(
        [
            facet_vertices[0]
            + sum((i - facet_vertices[0]) * j for i, j in zip(facet_vertices[1:], p))
            for p in points
        ],
        dtype=np.float64,
    )


def map_edge_points(
    points: npt.NDArray[np.float64], edge: int, cellname: str
) -> npt.NDArray[np.float64]:
    """Map points from a reference edge to a physical edge."""
    geom = np.asarray(basix.geometry(_CellType[cellname]))
    edge_vertices = [geom[i] for i in basix.topology(_CellType[cellname])[-3][edge]]
    return np.asarray(
        [
            edge_vertices[0] + sum((i - edge_vertices[0]) * j for i, j in zip(edge_vertices[1:], p))
            for p in points
        ],
        dtype=np.float64,
    )
