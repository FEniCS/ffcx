# Copyright (C) 2012-2017 Marie Rognes
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Utility functions for some code shared between representations."""

import hashlib
import logging

import numpy

import ufl
from ffcx.fiatinterface import (create_quadrature, map_facet_points,
                                reference_cell_vertices)

logger = logging.getLogger(__name__)


class QuadratureRule:
    def __init__(self, points, weights):
        self.points = points
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
        (points, weights) = create_quadrature(cell.cellname(), degree, rule)
    elif integral_type in ufl.measure.facet_integral_types:
        (points, weights) = create_quadrature(ufl.cell.cellname2facetname[cell.cellname()], degree,
                                              rule)
    elif integral_type in ufl.measure.point_integral_types:
        (points, weights) = create_quadrature("vertex", degree, rule)
    elif integral_type in ufl.custom_integral_types:
        (points, weights) = (None, None)
    elif integral_type == "expression":
        (points, weights) = (None, None)
    else:
        logging.exception("Unknown integral type: {}".format(integral_type))
    return (points, weights)


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
        raise RuntimeError("Unknown integral_type: {}".format(integral_type))
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
        raise RuntimeError("Can't map points from entity_dim=%s" % (entity_dim, ))


def generate_enabled_coefficients(enabled_coefficients):
    # TODO: I don't know how to implement this using the format dict,
    # this will do for now:
    initializer_list = ", ".join("true" if enabled else "false" for enabled in enabled_coefficients)
    if enabled_coefficients:
        code = '\n'.join(["[{}] = {{ {} }};".format(len(enabled_coefficients), initializer_list)])
    else:
        code = "[1] = {false};  /* No coefficients, but C does not permit zero-sized arrays */"
    return code


def initialize_integral_code(ir, parameters):
    """Default integral IR.

    Representation independent default initialization of code dict for
    integral from intermediate representation.

    """
    code = {}
    code["class_type"] = ir.integral_type + "_integral"
    code["name"] = ir.name
    code["members"] = ""
    code["constructor"] = ""
    code["constructor_arguments"] = ""
    code["initializer_list"] = ""
    code["destructor"] = ""
    code["enabled_coefficients"] = generate_enabled_coefficients(ir.enabled_coefficients)
    code["additional_includes_set"] = set()  # FIXME: Get this out of code[]

    return code


def initialize_expression_code(ir):
    code = {}
    code["name"] = "{}_expression".format(ir.name)

    return code
