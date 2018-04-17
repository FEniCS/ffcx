# -*- coding: utf-8 -*-
"""This module contains utility functions for some code shared between
representations.

"""

# Copyright (C) 2012-2017 Marie Rognes
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Martin Sandve Alnæs 2013-2017
# Modified by Anders Logg 2014

import logging

import numpy

from ffc.classname import make_integral_classname
from ffc.fiatinterface import (create_element, create_quadrature,
                               map_facet_points, reference_cell_vertices)
from ufl.cell import cellname2facetname
from ufl.measure import (custom_integral_types, facet_integral_types,
                         integral_type_to_measure_name, point_integral_types)

logger = logging.getLogger(__name__)


def create_quadrature_points_and_weights(integral_type, cell, degree, rule):
    """Create quadrature rule and return points and weights."""
    if integral_type == "cell":
        (points, weights) = create_quadrature(cell.cellname(), degree, rule)
    elif integral_type in facet_integral_types:
        (points, weights) = create_quadrature(cellname2facetname[cell.cellname()], degree, rule)
    elif integral_type in point_integral_types:
        (points, weights) = create_quadrature("vertex", degree, rule)
    elif integral_type in custom_integral_types:
        (points, weights) = (None, None)
    else:
        logging.exception("Unknown integral type: {}".format(integral_type))
    return (points, weights)


def integral_type_to_entity_dim(integral_type, tdim):
    """Given integral_type and domain tdim, return the tdim of the
    integration entity.

    """
    if integral_type == "cell":
        entity_dim = tdim
    elif integral_type in facet_integral_types:
        entity_dim = tdim - 1
    elif integral_type in point_integral_types:
        entity_dim = 0
    elif integral_type in custom_integral_types:
        entity_dim = tdim
    else:
        logger.exception("Unknown integral_type: {}".format(integral_type))
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
        logger.exception("Can't map points from entity_dim=%s" % (entity_dim, ))


def needs_oriented_jacobian(form_data):
    # Check whether this form needs an oriented jacobian (only forms
    # involgin contravariant piola mappings seem to need it)
    for ufl_element in form_data.unique_elements:
        element = create_element(ufl_element)
        if "contravariant piola" in element.mapping():
            return True
    return False


# Mapping from recognized domain types to entity types
_entity_types = {
    "cell": "cell",
    "exterior_facet": "facet",
    "interior_facet": "facet",
    "vertex": "vertex",
    "custom": "cell"
}


def entity_type_from_integral_type(integral_type):
    return _entity_types[integral_type]


def initialize_integral_ir(representation, itg_data, form_data, form_id):
    """Initialize a representation dict with common information that is
    expected independently of which representation is chosen.

    """

    entitytype = entity_type_from_integral_type(itg_data.integral_type)
    cell = itg_data.domain.ufl_cell()
    #cellname = cell.cellname()
    tdim = cell.topological_dimension()
    assert all(tdim == itg.ufl_domain().topological_dimension() for itg in itg_data.integrals)

    # Set number of cells if not set TODO: Get automatically from
    # number of domains
    num_cells = itg_data.metadata.get("num_cells")

    return {
        "representation": representation,
        "integral_type": itg_data.integral_type,
        "subdomain_id": itg_data.subdomain_id,
        "form_id": form_id,
        "rank": form_data.rank,
        "geometric_dimension": form_data.geometric_dimension,
        "topological_dimension": tdim,
        "entitytype": entitytype,
        "num_facets": cell.num_facets(),
        "num_vertices": cell.num_vertices(),
        "needs_oriented": needs_oriented_jacobian(form_data),
        "num_cells": num_cells,
        "enabled_coefficients": itg_data.enabled_coefficients
    }


def generate_enabled_coefficients(enabled_coefficients):
    # TODO: I don't know how to implement this using the format dict,
    # this will do for now:
    initializer_list = ", ".join("true" if enabled else "false" for enabled in enabled_coefficients)
    if enabled_coefficients:
        code = '\n'.join(["[{}] = {{ {} }};".format(len(enabled_coefficients), initializer_list)])
    else:
        code = "[] = {};"
    return code


def initialize_integral_code(ir, prefix, parameters):
    """Representation independent default initialization of code dict for
    integral from intermediate representation.

    """
    code = {}
    code["class_type"] = ir["integral_type"] + "_integral"
    code["classname"] = make_integral_classname(prefix, ir["integral_type"], ir["form_id"],
                                                ir["subdomain_id"])
    code["members"] = ""
    code["constructor"] = ""
    code["constructor_arguments"] = ""
    code["initializer_list"] = ""
    code["destructor"] = ""
    code["enabled_coefficients"] = generate_enabled_coefficients(ir["enabled_coefficients"])
    code["additional_includes_set"] = set()  # FIXME: Get this out of code[]

    return code
