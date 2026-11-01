# Copyright (C) 2021 Matthew W. Scroggs and Chris Richardson
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Finite element interface."""

import typing

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


def tensor_product_factorisation(
    element: basix.ufl._ElementBase,
) -> tuple[list[basix.finite_element.FiniteElement], npt.NDArray[np.int32] | None] | None:
    """Get the 1D factors of an element's basis and the dofs they give.

    The factors are the elements whose basis functions, multiplied
    together with the first direction's dof index varying slowest, give
    the element's basis. Returns the factors and the element dof of each
    such product, or None for the dofs of an element whose own dofs are
    numbered this way. Returns None if the element does not factorise.
    """
    if isinstance(element, basix.ufl._BlockedElement):
        element = typing.cast(basix.ufl._ElementBase, element.sub_elements[0])
    if element.has_tensor_product_factorisation:
        factors = element.get_tensor_product_representation()
        return (factors[0], None) if factors is not None and len(factors) == 1 else None
    if isinstance(element, basix.ufl._BasixElement):
        return _ordered_tensor_product_factorisation(element)
    return None


_ordered_factorisations: dict[basix.ufl._BasixElement, tuple | None] = {}


def _ordered_tensor_product_factorisation(
    element: basix.ufl._BasixElement,
) -> tuple[list[basix.finite_element.FiniteElement], npt.NDArray[np.int32]] | None:
    """Factorise an element numbering its dofs by entity, such as Lagrange.

    The element with tensor product ordering of the same family and
    variants supplies the factors. The dof of each product of factors is
    found by matching basis functions at points in the cell.
    """
    if element in _ordered_factorisations:
        return _ordered_factorisations[element]
    e = element.basix_element
    result = None
    if (
        e.cell_type in (_CellType.quadrilateral, _CellType.hexahedron)
        and e.value_size == 1
        and len(e.dof_ordering) == 0
    ):
        try:
            tp = basix.create_tp_element(
                e.family, e.cell_type, e.degree, e.lagrange_variant, e.dpc_variant, e.discontinuous
            )
            factors = (
                tp.get_tensor_product_representation()
                if tp.has_tensor_product_factorisation
                else []
            )
        except RuntimeError:
            factors = []
        if len(factors) == 1:
            dim = len(factors[0])
            points = np.random.default_rng(0).random((4 * e.dim, dim))
            values = e.tabulate(0, points)[0, :, :, 0]
            product = np.ones((points.shape[0], 1))
            for d, factor in enumerate(factors[0]):
                f = factor.tabulate(0, points[:, d : d + 1])[0, :, :, 0]
                product = (product[:, :, None] * f[:, None, :]).reshape(points.shape[0], -1)
            if product.shape == values.shape:
                scale = np.abs(values).max()
                match = np.isclose(
                    product[:, :, None], values[:, None, :], rtol=0, atol=1e-10 * scale
                ).all(axis=0)
                dofs = np.flatnonzero(match.any(axis=1))
                if len(dofs) == product.shape[1] and (match.sum(axis=1) == 1).all():
                    perm = match.argmax(axis=1).astype(np.int32)
                    if len(np.unique(perm)) == len(perm):
                        result = (factors[0], perm)
    _ordered_factorisations[element] = result
    return result


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
