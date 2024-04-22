# Copyright (C) 2012-2017 Marie Rognes
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

from ffcx.element_interface import create_quadrature, map_facet_points, reference_cell_vertices

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
    integral_type, cell, degree, rule, elements, use_tensor_product=False
):
    """Create quadrature rule and return points and weights."""
    pts = None
    wts = None
    tensor_factors = None

    if integral_type == "cell":
        if cell.cellname() in ["quadrilateral", "hexahedron"] and use_tensor_product:
            if cell.cellname() == "quadrilateral":
                tensor_factors = [
                    create_quadrature("interval", degree, rule, elements) for _ in range(2)
                ]
            elif cell.cellname() == "hexahedron":
                tensor_factors = [
                    create_quadrature("interval", degree, rule, elements) for _ in range(3)
                ]
            pts = np.array(
                [tuple(i[0] for i in p) for p in itertools.product(*[f[0] for f in tensor_factors])]
            )
            wts = np.array([np.prod(p) for p in itertools.product(*[f[1] for f in tensor_factors])])
        else:
            pts, wts = create_quadrature(cell.cellname(), degree, rule, elements)
    elif integral_type in ufl.measure.facet_integral_types:
        facet_types = cell.facet_types()
        # Raise exception for cells with more than one facet type e.g. prisms
        if len(facet_types) > 1:
            raise Exception(f"Cell type {cell} not supported for integral type {integral_type}.")
        pts, wts = create_quadrature(facet_types[0].cellname(), degree, rule, elements)
    elif integral_type in ufl.measure.point_integral_types:
        pts, wts = create_quadrature("vertex", degree, rule, elements)
    elif integral_type == "expression":
        pass
    else:
        logging.exception(f"Unknown integral type: {integral_type}")

    return pts, wts, tensor_factors


def integral_type_to_entity_dim(integral_type, tdim):
    """Given integral_type and domain tdim, return the tdim of the integration entity."""
    if integral_type == "cell":
        entity_dim = tdim
    elif integral_type in ufl.measure.facet_integral_types:
        entity_dim = tdim - 1
    elif integral_type in ufl.measure.point_integral_types:
        entity_dim = 0
    elif integral_type in ufl.custom_integral_types:
        entity_dim = tdim
    elif integral_type == "expression":
        entity_dim = tdim
    else:
        raise RuntimeError(f"Unknown integral_type: {integral_type}")
    return entity_dim


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
    elif entity_dim == 0:
        return np.asarray([reference_cell_vertices(cell.cellname())[entity]])
    else:
        raise RuntimeError(f"Can't map points from entity_dim={entity_dim}")
