# -*- coding: utf-8 -*-
"""This module contains utility functions for some code shared between
quadrature and tensor representation."""

# Copyright (C) 2012-2015 Marie Rognes
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
# Modified by Martin Alnaes 2013-2015
# Modified by Anders Logg 2014

import numpy

from ufl.measure import integral_type_to_measure_name
from ufl.cell import cellname2facetname
from ufl import custom_integral_types

from ffc.fiatinterface import create_element
from ffc.cpp import make_integral_classname
from ffc.log import error

from ffc.fiatinterface import create_quadrature


def create_quadrature_points_and_weights(integral_type, cell, degree, rule):
    "Create quadrature rule and return points and weights."
    if integral_type == "cell":
        (points, weights) = create_quadrature(cell.cellname(), degree, rule)
    elif integral_type == "exterior_facet" or integral_type == "interior_facet":
        facet_cellname = cellname2facetname[cell.cellname()]
        (points, weights) = create_quadrature(facet_cellname, degree, rule)
    elif integral_type == "vertex":
        (points, weights) = ([()], numpy.array([1.0, ]))  # TODO: Will be fixed
    elif integral_type in custom_integral_types:
        (points, weights) = (None, None)
    else:
        error("Unknown integral type: " + str(integral_type))
    return (points, weights)


def transform_component(component, offset, ufl_element):
    """
    This function accounts for the fact that if the geometrical and
    topological dimension does not match, then for native vector
    elements, in particular the Piola-mapped ones, the physical value
    dimensions and the reference value dimensions are not the
    same. This has certain consequences for mixed elements, aka 'fun
    with offsets'.
    """
    # This code is used for tensor/monomialtransformation.py and
    # quadrature/quadraturetransformerbase.py.

    cell = ufl_element.cell()
    gdim = cell.geometric_dimension()
    tdim = cell.topological_dimension()

    # Do nothing if we are not in a special case: The special cases
    # occur if we have piola mapped elements (for which value_shape !=
    # ()), and if gdim != tdim)
    if gdim == tdim:
        return component, offset
    all_mappings = create_element(ufl_element).mapping()
    special_case = (any(['piola' in m for m in all_mappings])
                    and ufl_element.num_sub_elements() > 1)
    if not special_case:
        return component, offset

    # Extract lists of reference and physical value dimensions by
    # sub-element
    reference_value_dims = []
    physical_value_dims = []
    for sub_element in ufl_element.sub_elements():
        assert (len(sub_element.value_shape()) < 2), \
            "Vector-valued assumption failed"
        if sub_element.value_shape() == ():
            reference_value_dims += [1]
            physical_value_dims += [1]
        else:
            reference_value_dims += [sub_element.value_shape()[0]
                                     - (gdim - tdim)]
            physical_value_dims += [sub_element.value_shape()[0]]

    # Figure out which sub-element number 'component' is in,
    # 'sub_element_number' contains the result
    tot = physical_value_dims[0]
    for sub_element_number in range(len(physical_value_dims)):
        if component < tot:
            break
        else:
            tot += physical_value_dims[sub_element_number + 1]

    # Compute the new reference offset:
    reference_offset = sum(reference_value_dims[:sub_element_number])
    physical_offset = sum(physical_value_dims[:sub_element_number])
    shift = physical_offset - reference_offset

    # Compute the component relative to the reference frame
    reference_component = component - shift

    return reference_component, reference_offset


def needs_oriented_jacobian(form_data):
    # Check whether this form needs an oriented jacobian (only forms
    # involgin contravariant piola mappings seem to need it)
    for ufl_element in form_data.unique_elements:
        element = create_element(ufl_element)
        if "contravariant piola" in element.mapping():
            return True
    return False


def initialize_integral_ir(representation, itg_data, form_data, form_id):
    """Initialize a representation dict with common information that is
    expected independently of which representation is chosen."""

    # Mapping from recognized domain types to entity types
    entity_type = {"cell": "cell",
                   "exterior_facet": "facet",
                   "interior_facet": "facet",
                   "vertex": "vertex",
                   # "point":          "vertex", # TODO: Not sure, clarify here what 'entity_type' refers to?
                   "custom": "cell",
                   "cutcell": "cell",
                   "interface": "cell",
                   "overlap": "cell",
                   }[itg_data.integral_type]

    # Extract data
    cell = itg_data.domain.ufl_cell()
    cellname = cell.cellname()
    tdim = cell.topological_dimension()
    assert all(tdim == itg.ufl_domain().topological_dimension() for itg in itg_data.integrals)

    # Set number of cells if not set TODO: Get automatically from number of domains
    num_cells = itg_data.metadata.get("num_cells")

    return {"representation": representation,
            "integral_type": itg_data.integral_type,
            "subdomain_id": itg_data.subdomain_id,
            "form_id": form_id,
            "rank": form_data.rank,
            "geometric_dimension": form_data.geometric_dimension,
            "topological_dimension": tdim,
            "entitytype": entity_type,
            "num_facets": cell.num_facets(),
            "num_vertices": cell.num_vertices(),
            "needs_oriented": needs_oriented_jacobian(form_data),
            "num_cells": num_cells,
            "enabled_coefficients": itg_data.enabled_coefficients,
            }


def generate_enabled_coefficients(enabled_coefficients):
    # TODO: I don't know how to implement this using the format dict, this will do for now:
    initializer_list = ", ".join("true" if enabled else "false"
                                 for enabled in enabled_coefficients)
    code = '\n'.join([
        "static const std::vector<bool> enabled({%s});" % initializer_list,
        "return enabled;",
    ])
    return code


def initialize_integral_code(ir, prefix, parameters):
    "Representation independent default initialization of code dict for integral from intermediate representation."
    code = {}
    code["class_type"] = ir["integral_type"] + "_integral"
    code["classname"] = make_integral_classname(prefix, ir["integral_type"], ir["form_id"], ir["subdomain_id"])
    code["members"] = ""
    code["constructor"] = ""
    code["constructor_arguments"] = ""
    code["initializer_list"] = ""
    code["destructor"] = ""
    code["enabled_coefficients"] = generate_enabled_coefficients(ir["enabled_coefficients"])
    # code["additional_includes_set"] = set() #ir["additional_includes_set"]
    return code
