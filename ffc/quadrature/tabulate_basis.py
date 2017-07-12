# -*- coding: utf-8 -*-

# Copyright (C) 2009-2014 Kristian B. Oelgaard
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
# Modified by Anders Logg, 2009, 2015
# Modified by Martin Sandve Alnæs, 2013-2014

"Quadrature representation class."

import numpy
import itertools

# UFL modules
from ufl.cell import num_cell_entities
from ufl.classes import ReferenceGrad, Grad, CellAvg, FacetAvg
from ufl.algorithms import extract_unique_elements, extract_type, extract_elements
from ufl import custom_integral_types

# FFC modules
from ffc.log import error
from ffc.utils import product, insert_nested_dict
from ffc.fiatinterface import create_element
from ffc.fiatinterface import map_facet_points, reference_cell_vertices
from ffc.representationutils import create_quadrature_points_and_weights
from ffc.representationutils import integral_type_to_entity_dim


def _find_element_derivatives(expr, elements, element_replace_map):
    "Find the highest derivatives of given elements in expression."
    # TODO: This is most likely not the best way to get the highest
    #       derivative of an element, but it works!

    # Initialise dictionary of elements and the number of derivatives.
    # (Note that elements are already mapped through the
    # element_replace_map)
    num_derivatives = dict((e, 0) for e in elements)

    # Extract the derivatives from the integral.
    derivatives = set(extract_type(expr, Grad)) | set(extract_type(expr,
                                                                   ReferenceGrad))

    # Loop derivatives and extract multiple derivatives.
    for d in list(derivatives):
        # After UFL has evaluated derivatives, only one element can be
        # found inside any single Grad expression
        elem, = extract_elements(d.ufl_operands[0])
        elem = element_replace_map[elem]
        # Set the number of derivatives to the highest value
        # encountered so far.
        num_derivatives[elem] = max(num_derivatives[elem],
                                    len(extract_type(d, Grad)),
                                    len(extract_type(d, ReferenceGrad)))
    return num_derivatives


def _tabulate_empty_psi_table(tdim, deriv_order, element):
    "Tabulate psi table when there are no points (custom integrals)."

    # All combinations of partial derivatives up to given order
    gdim = tdim  # hack, consider passing gdim variable here
    derivs = [d for d in itertools.product(*(gdim*[list(range(0, deriv_order + 1))]))]
    derivs = [d for d in derivs if sum(d) <= deriv_order]

    # Return empty table
    table = {}
    for d in derivs:
        value_shape = element.value_shape()
        if value_shape == ():
            table[d] = [[]]
        else:
            value_size = product(value_shape)
            table[d] = [[[] for c in range(value_size)]]

    # Let entity be 0 even for non-cells, this is for
    # custom integrals where we don't need tables to
    # contain multiple entitites
    entity = 0
    return {entity: table}


def _map_entity_points(cellname, tdim, points, entity_dim, entity):
    # Not sure if this is useful anywhere else than in _tabulate_psi_table!
    if entity_dim == tdim:
        assert entity == 0
        return points
    elif entity_dim == tdim-1:
        return map_facet_points(points, entity, cellname)
    elif entity_dim == 0:
        return (reference_cell_vertices(cellname)[entity],)


def _tabulate_psi_table(integral_type, cellname, tdim,
                        element, deriv_order, points):
    "Tabulate psi table for different integral types."
    # Handle case when list of points is empty
    if points is None:
        return _tabulate_empty_psi_table(tdim, deriv_order, element)

    # Otherwise, call FIAT to tabulate
    entity_dim = integral_type_to_entity_dim(integral_type, tdim)
    num_entities = num_cell_entities[cellname][entity_dim]
    psi_table = {}
    for entity in range(num_entities):
        entity_points = _map_entity_points(cellname, tdim, points, entity_dim, entity)
        psi_table[entity] = element.tabulate(deriv_order, entity_points)
    return psi_table


# MSA: This function is in serious need for some refactoring and
#      splitting up.  Or perhaps I should just add a new
#      implementation for uflacs, but I'd rather not have two versions
#      to maintain.
def tabulate_basis(sorted_integrals, form_data, itg_data):
    "Tabulate the basisfunctions and derivatives."

    # MER: Note to newbies: this code assumes that each integral in
    # the dictionary of sorted_integrals that enters here, has a
    # unique number of quadrature points ...

    # Initialise return values.
    quadrature_rules = {}
    psi_tables = {}
    integrals = {}
    avg_elements = {"cell": [], "facet": []}

    # Get some useful variables in short form
    integral_type = itg_data.integral_type
    cell = itg_data.domain.ufl_cell()
    cellname = cell.cellname()
    tdim = itg_data.domain.topological_dimension()
    entity_dim = integral_type_to_entity_dim(integral_type, tdim)
    num_entities = num_cell_entities[cellname][entity_dim]

    # Create canonical ordering of quadrature rules
    rules = sorted(sorted_integrals.keys())

    # Loop the quadrature points and tabulate the basis values.
    for rule in rules:
        scheme, degree = rule

        # --------- Creating quadrature rule
        # Make quadrature rule and get points and weights.
        (points, weights) = create_quadrature_points_and_weights(integral_type,
                                                                 cell, degree,
                                                                 scheme)

        # The TOTAL number of weights/points
        num_points = None if weights is None else len(weights)

        # Add points and rules to dictionary
        if num_points in quadrature_rules:
            error("This number of points is already present in the weight table: " + repr(quadrature_rules))
        quadrature_rules[num_points] = (weights, points)

        # --------- Store integral
        # Add the integral with the number of points as a key to the
        # return integrals.
        integral = sorted_integrals[rule]
        if num_points in integrals:
            error("This number of points is already present in the integrals: " + repr(integrals))
        integrals[num_points] = integral

        # --------- Analyse UFL elements in integral

        # Get all unique elements in integral.
        ufl_elements = [form_data.element_replace_map[e]
                        for e in extract_unique_elements(integral)]

        # Insert elements for x and J
        domain = integral.ufl_domain()  # FIXME: For all domains to be sure? Better to rewrite though.
        x_element = domain.ufl_coordinate_element()
        if x_element not in ufl_elements:
            if integral_type in custom_integral_types:
                # FIXME: Not yet implemented, in progress
                # warning("Vector elements not yet supported in custom integrals so element for coordinate function x will not be generated.")
                pass
            else:
                ufl_elements.append(x_element)

        # Find all CellAvg and FacetAvg in integrals and extract
        # elements
        for avg, AvgType in (("cell", CellAvg), ("facet", FacetAvg)):
            expressions = extract_type(integral, AvgType)
            avg_elements[avg] = [form_data.element_replace_map[e]
                                 for expr in expressions
                                 for e in extract_unique_elements(expr)]

        # Find the highest number of derivatives needed for each element
        num_derivatives = _find_element_derivatives(integral.integrand(),
                                                    ufl_elements,
                                                    form_data.element_replace_map)
        # Need at least 1 for the Jacobian
        num_derivatives[x_element] = max(num_derivatives.get(x_element, 0), 1)

        # --------- Evaluate FIAT elements in quadrature points and
        # --------- store in tables

        # Add the number of points to the psi tables dictionary
        if num_points in psi_tables:
            error("This number of points is already present in the psi table: " + repr(psi_tables))
        psi_tables[num_points] = {}

        # Loop FIAT elements and tabulate basis as usual.
        for ufl_element in ufl_elements:
            fiat_element = create_element(ufl_element)

            # Tabulate table of basis functions and derivatives in
            # points
            psi_table = _tabulate_psi_table(integral_type, cellname, tdim,
                                            fiat_element,
                                            num_derivatives[ufl_element],
                                            points)

            # Insert table into dictionary based on UFL elements
            # (None=not averaged)
            avg = None
            psi_tables[num_points][ufl_element] = { avg: psi_table }

    # Loop over elements found in CellAvg and tabulate basis averages
    num_points = 1
    for avg in ("cell", "facet"):
        # Doesn't matter if it's exterior or interior
        if avg == "cell":
            avg_integral_type = "cell"
        elif avg == "facet":
            avg_integral_type = "exterior_facet"

        for element in avg_elements[avg]:
            fiat_element = create_element(element)

            # Make quadrature rule and get points and weights.
            (points, weights) = create_quadrature_points_and_weights(avg_integral_type, cell, element.degree(), "default")
            wsum = sum(weights)

            # Tabulate table of basis functions and derivatives in
            # points
            entity_psi_tables = _tabulate_psi_table(avg_integral_type,
                                                    cellname, tdim,
                                                    fiat_element, 0, points)
            rank = len(element.value_shape())

            # Hack, duplicating table with per-cell values for each
            # facet in the case of cell_avg(f) in a facet integral
            if num_entities > len(entity_psi_tables):
                assert len(entity_psi_tables) == 1
                assert avg_integral_type == "cell"
                assert "facet" in integral_type
                v, = sorted(entity_psi_tables.values())
                entity_psi_tables = dict((e, v) for e in range(num_entities))

            for entity, deriv_table in sorted(entity_psi_tables.items()):
                deriv, = sorted(deriv_table.keys())  # Not expecting derivatives of averages
                psi_table = deriv_table[deriv]

                if rank:
                    # Compute numeric integral
                    num_dofs, num_components, num_points = psi_table.shape
                    if num_points != len(weights):
                        error("Weights and table shape does not match.")
                    avg_psi_table = numpy.asarray([[[numpy.dot(psi_table[j, k, :], weights) / wsum]
                                                    for k in range(num_components)]
                                                   for j in range(num_dofs)])
                else:
                    # Compute numeric integral
                    num_dofs, num_points = psi_table.shape
                    if num_points != len(weights):
                        error("Weights and table shape does not match.")
                    avg_psi_table = numpy.asarray([[numpy.dot(psi_table[j, :],
                                                              weights) / wsum] for j in range(num_dofs)])

                # Insert table into dictionary based on UFL elements
                insert_nested_dict(psi_tables, (num_points, element, avg,
                                                entity, deriv), avg_psi_table)

    return (integrals, psi_tables, quadrature_rules)
