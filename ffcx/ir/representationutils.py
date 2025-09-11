# Copyright (C) 2012-2025 Marie Rognes and Jørgen S. Dokken
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Utility functions for some code shared between representations."""

import hashlib
import itertools
import logging

import numpy as np
import ufl
from basix.ufl import _ElementBase

from ffcx.definitions import IntegralType
from ffcx.element_interface import (
    create_quadrature,
    map_edge_points,
    map_facet_points,
    reference_cell_vertices,
)

logger = logging.getLogger("ffcx")


class QuadratureRule:
    """A quadrature rule."""

    def __init__(self, points, weights, tensor_factors=None):
        """Initialise."""
        self.points = np.ascontiguousarray(points)  # TODO: change basix to make this unnecessary
        self.weights = weights
        self.tensor_factors = tensor_factors
        self.has_tensor_factors = tensor_factors is not None
        self._hash = None

    def __hash__(self):
        """Hash."""
        if self._hash is None:
            self.hash_obj = hashlib.sha1(self.points)
            self._hash = int(self.hash_obj.hexdigest(), 32)
        return self._hash

    def __eq__(self, other):
        """Check equality."""
        return np.allclose(self.points, other.points) and np.allclose(self.weights, other.weights)

    def id(self):
        """Return unique deterministic identifier.

        Note:
            This identifier is used to provide unique names to tables and symbols
            in generated code.
        """
        return self.hash_obj.hexdigest()[-3:]


def create_quadrature_points_and_weights(
    integral_type: IntegralType,
    cell: ufl.Cell,
    degree: int,
    rule: str,
    elements: list[_ElementBase],
    use_tensor_product: bool = False,
):
    """Create quadrature rule and return points and weights."""
    pts = {}
    wts = {}
    tensor_factors = {}
    match integral_type:
        case IntegralType(is_expression=True):
            pass
        case IntegralType(codim=0):
            cell_name = cell.cellname()
            if cell_name in ["quadrilateral", "hexahedron"] and use_tensor_product:
                if cell_name == "quadrilateral":
                    tensor_factors[cell_name] = [
                        create_quadrature("interval", degree, rule, elements) for _ in range(2)
                    ]
                elif cell_name == "hexahedron":
                    tensor_factors[cell_name] = [
                        create_quadrature("interval", degree, rule, elements) for _ in range(3)
                    ]
                pts[cell_name] = np.array(
                    [
                        tuple(i[0] for i in p)  # type: ignore
                        for p in itertools.product(*[f[0] for f in tensor_factors[cell_name]])  # type: ignore
                    ]
                )
                wts[cell_name] = np.array(
                    [
                        np.prod(p)  # type: ignore
                        for p in itertools.product(*[f[1] for f in tensor_factors[cell_name]])  # type: ignore
                    ]
                )
            else:
                pts[cell_name], wts[cell_name] = create_quadrature(  # type: ignore
                    cell_name, degree, rule, elements
                )
        case IntegralType(codim=1):
            for ft in cell.facet_types():
                pts[ft.cellname()], wts[ft.cellname()] = create_quadrature(  # type: ignore
                    ft.cellname(),
                    degree,
                    rule,
                    elements,
                )
        case IntegralType(codim=2):
            for rt in cell.ridge_types():
                pts[rt.cellname()], wts[rt.cellname()] = create_quadrature(  # type: ignore
                    rt.cellname(),
                    degree,
                    rule,
                    elements,
                )
        case IntegralType(codim=-1):
            pts["vertex"], wts["vertex"] = create_quadrature("vertex", degree, rule, elements)  # type: ignore
        case _:
            raise RuntimeError(
                f"Unknown integral_type {integral_type} in quadrature rule creation."
            )
    return pts, wts, tensor_factors


def integral_type_to_entity_dim(integral_type: IntegralType, tdim: int):
    """Given integral_type and domain tdim, return the tdim of the integration entity."""
    if integral_type.codim == -1:
        return 0
    elif integral_type.codim >= 0:
        return tdim - integral_type.codim
    else:
        raise NotImplementedError(f"Cannot handle codim {integral_type.codim} < -1")


def map_integral_points(points, integral_type, cell, entity):
    """Map points from reference entity to its parent reference cell."""
    tdim = cell.topological_dimension()
    entity_dim = integral_type_to_entity_dim(integral_type, tdim)
    if entity_dim == tdim:
        assert points.shape[1] == tdim
        assert entity == 0
        return np.asarray(points)
    elif entity_dim == tdim - 1:
        assert points.shape[1] == tdim - 1
        return np.asarray(map_facet_points(points, entity, cell.cellname()))
    elif entity_dim == tdim - 2:
        assert points.shape[1] == tdim - 2
        # Special handling of pushing forward 0D points to cell
        if entity_dim == 0:
            assert points.shape[1] == 0
            points = np.zeros((1, 1))
        return np.asarray(map_edge_points(points, entity, cell.cellname()))
    elif entity_dim == 0:
        return np.asarray([reference_cell_vertices(cell.cellname())[entity]])
    else:
        raise RuntimeError(f"Can't map points from entity_dim={entity_dim}")
