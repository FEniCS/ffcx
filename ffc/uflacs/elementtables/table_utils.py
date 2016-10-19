# -*- coding: utf-8 -*-
# Copyright (C) 2011-2016 Martin Sandve Alnæs
#
# This file is part of UFLACS.
#
# UFLACS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# UFLACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with UFLACS. If not, see <http://www.gnu.org/licenses/>

"""Utilities for precomputed table manipulation."""

from __future__ import print_function  # used in some debugging

import numpy

from ufl.permutation import build_component_numbering
from ufl.cell import num_cell_entities

from ffc.log import error
from ffc.fiatinterface import create_element
from ffc.representationutils import integral_type_to_entity_dim, map_integral_points
from ffc.representationutils import create_quadrature_points_and_weights

def equal_tables(a, b, eps):
    "Compare tables to be equal within a tolerance."
    a = numpy.asarray(a)
    b = numpy.asarray(b)
    if a.shape != b.shape:
        return False
    if len(a.shape) > 1:
        return all(equal_tables(a[i], b[i], eps)
                   for i in range(a.shape[0]))
    def scalars_equal(x, y, eps):
        return abs(x-y) < eps
    return all(scalars_equal(a[i], b[i], eps)
               for i in range(a.shape[0]))


def clamp_table_small_integers(table, eps):
    "Clamp almost 0,1,-1 values to integers. Returns new table."
    # Get shape of table and number of columns, defined as the last axis
    table = numpy.asarray(table)
    for n in (-1, 0, 1):
        table[numpy.where(abs(table - n) < eps)] = float(n)
    return table


def strip_table_zeros(table, eps):
    "Strip zero columns from table. Returns column range (begin,end) and the new compact table."
    # Get shape of table and number of columns, defined as the last axis
    table = numpy.asarray(table)
    sh = table.shape
    nc = sh[-1]

    # Find first nonzero column
    begin = nc
    for i in range(nc):
        if numpy.linalg.norm(table[..., i]) > eps:
            begin = i
            break

    # Find (one beyond) last nonzero column
    end = begin
    for i in range(nc-1, begin-1, -1):
        if numpy.linalg.norm(table[..., i]) > eps:
            end = i+1
            break

    # Make subtable by stripping first and last columns
    stripped_table = table[..., begin:end]
    return begin, end, stripped_table


def build_unique_tables(tables, eps):
    """Given a list or dict of tables, return a list of unique tables
    and a dict of unique table indices for each input table key."""
    unique = []
    mapping = {}

    if isinstance(tables, list):
        keys = list(range(len(tables)))
    elif isinstance(tables, dict):
        keys = sorted(tables.keys())

    for k in keys:
        t = tables[k]
        found = -1
        for i, u in enumerate(unique):
            if equal_tables(u, t, eps):
                found = i
                break
        if found == -1:
            i = len(unique)
            unique.append(t)
        mapping[k] = i

    return unique, mapping


def get_ffc_table_values(points,
                         cell, integral_type,
                         num_points, # TODO: Remove, not needed
                         ufl_element, avg,
                         entitytype, derivative_counts,
                         flat_component, epsilon):
    """Extract values from ffc element table.

    Returns a 3D numpy array with axes
    (entity number, quadrature point number, dof number)
    """
    deriv_order = sum(derivative_counts)

    if avg in ("cell", "facet"):
        # Redefine points to compute average tables

        # Make sure this is not called with points, that doesn't make sense
        #assert points is None
        #assert num_points is None

        # Not expecting derivatives of averages
        assert not any(derivative_counts)
        assert deriv_order == 0

        # Doesn't matter if it's exterior or interior facet integral,
        # just need a valid integral type to create quadrature rule
        if avg == "cell":
            integral_type = "cell"
        elif avg == "facet":
            integral_type = "exterior_facet"

        # Make quadrature rule and get points and weights
        points, weights = create_quadrature_points_and_weights(
            integral_type, cell, ufl_element.degree(), "default")

    # Tabulate table of basis functions and derivatives in points for each entity
    fiat_element = create_element(ufl_element)
    tdim = cell.topological_dimension()
    entity_dim = integral_type_to_entity_dim(integral_type, tdim)
    num_entities = num_cell_entities[cell.cellname()][entity_dim]
    entity_tables = []
    for entity in range(num_entities):
        entity_points = map_integral_points(points, integral_type, cell, entity)
        tbl = fiat_element.tabulate(deriv_order, entity_points)[derivative_counts]
        entity_tables.append(tbl)

    # Extract arrays for the right scalar component
    component_tables = []
    sh = ufl_element.value_shape()
    if sh == ():
        # Scalar valued element
        for entity, entity_table in enumerate(entity_tables):
            component_tables.append(entity_table)
    elif len(sh) == 2 and ufl_element.num_sub_elements() == 0:
        # 2-tensor-valued elements, not a tensor product
        # mapping flat_component back to tensor component
        (_, f2t) = build_component_numbering(sh, ufl_element.symmetry())
        t_comp = f2t[flat_component]
        for entity, entity_table in enumerate(entity_tables):
            tbl = entity_table[:, t_comp[0], t_comp[1], :]
            component_tables.append(tbl)
    else:
        # Vector-valued or mixed element
        for entity, entity_table in enumerate(entity_tables):
            tbl = entity_table[:, flat_component, :]
            component_tables.append(tbl)

    if avg in ("cell", "facet"):
        # Compute numeric integral of the each component table
        wsum = sum(weights)
        for entity, tbl in enumerate(component_tables):
            num_dofs = tbl.shape[0]
            tbl = numpy.dot(tbl, weights) / wsum
            tbl = numpy.reshape(tbl, (num_dofs, 1))
            component_tables[entity] = tbl

    # Loop over entities and fill table blockwise (each block = points x dofs)
    # Reorder axes as (points, dofs) instead of (dofs, points)
    assert len(component_tables) == num_entities
    num_dofs, num_points = component_tables[0].shape
    shape = (num_entities, num_points, num_dofs)
    res = numpy.zeros(shape)
    for entity in range(num_entities):
        res[entity, :, :] = numpy.transpose(component_tables[entity])
    return res


def generate_psi_table_name(num_points, element_counter, averaged,
                            entitytype, derivative_counts, flat_component):
    """Generate a name for the psi table of the form:
    FE#_C#_D###[_AC|_AF|][_F|V][_Q#], where '#' will be an integer value.

    FE  - is a simple counter to distinguish the various bases, it will be
          assigned in an arbitrary fashion.

    C   - is the component number if any (this does not yet take into account
          tensor valued functions)

    D   - is the number of derivatives in each spatial direction if any.
          If the element is defined in 3D, then D012 means d^3(*)/dydz^2.

    AC  - marks that the element values are averaged over the cell

    AF  - marks that the element values are averaged over the facet

    F   - marks that the first array dimension enumerates facets on the cell

    V   - marks that the first array dimension enumerates vertices on the cell

    Q   - number of quadrature points, to distinguish between tables in a mixed quadrature degree setting

    """
    name = "FE%d" % element_counter
    if flat_component is not None:
        name += "_C%d" % flat_component
    if any(derivative_counts):
        name += "_D" + "".join(str(d) for d in derivative_counts)
    name += { None: "", "cell": "_AC", "facet": "_AF" }[averaged]
    name += { "cell": "", "facet": "_F", "vertex": "_V" }[entitytype]
    if num_points is not None:
        name += "_Q%d" % num_points
    return name
