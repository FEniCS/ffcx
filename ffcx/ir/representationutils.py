# Copyright (C) 2012-2017 Marie Rognes
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Utility functions for some code shared between representations."""

import hashlib
import logging

import numpy
from ffcx.libtab_interface import create_quadrature, reference_cell_vertices, map_facet_points

import ufl

logger = logging.getLogger("ffcx")


class QuadratureRule:
    def __init__(self, points, weights):
        self.points = numpy.ascontiguousarray(points)  # TODO: change libtab to make this unnecessary
        self.weights = weights
        self._hash = None

    def __hash__(self):
        if self._hash is None:
            self.hash_obj = hashlib.sha1(self.points)
            self._hash = int(self.hash_obj.hexdigest(), 32)
        return self._hash

    def __eq__(self, other):
        return numpy.allclose(self.points, other.points) and numpy.allclose(self.weights, other.weights)

    def id(self):
        """Returns unique deterministic identifier.

        Note
        ----
        This identifier is used to provide unique names to tables and symbols
        in generated code.

        """
        return self.hash_obj.hexdigest()[-3:]


def create_quadrature_points_and_weights(integral_type, cell, degree, rule):
    """Create quadrature rule and return points and weights."""

    if integral_type == "cell":
        return create_quadrature(cell.cellname(), degree + 1, rule)
    elif integral_type in ufl.measure.facet_integral_types:
        return create_quadrature(ufl.cell.cellname2facetname[cell.cellname()], degree + 1, rule)
    elif integral_type in ufl.measure.point_integral_types:
        return create_quadrature("vertex", degree + 1, rule)
    elif integral_type == "expression":
        return (None, None)

    logging.exception(f"Unknown integral type: {integral_type}")
    return (None, None)


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
        return numpy.asarray(points)
    elif entity_dim == tdim - 1:
        assert points.shape[1] == tdim - 1
        return numpy.asarray(map_facet_points(points, entity, cell.cellname()))
    elif entity_dim == 0:
        return numpy.asarray([reference_cell_vertices(cell.cellname())[entity]])
    else:
        raise RuntimeError(f"Can't map points from entity_dim={entity_dim}")
