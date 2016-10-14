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

from six import itervalues, iterkeys
from six import advance_iterator as next

from ffc.log import error
from ffc.fiatinterface import create_element

from ufl.permutation import build_component_numbering


def equal_tables(a, b, eps):
    "Compare tables to be equal within a tolerance."
    a = numpy.asarray(a)
    b = numpy.asarray(b)
    if a.shape != b.shape:
        return False
    if len(a.shape) > 1:
        return all(equal_tables(a[i], b[i], eps) for i in range(a.shape[0]))

    def scalars_equal(x, y, eps):
        return abs(x-y) < eps
    return all(scalars_equal(a[i], b[i], eps) for i in range(a.shape[0]))


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

'''
# Straight copy from quadrature/tabulate_basis
def _tabulate_psi_table(integral_type, cellname, tdim, fiat_element, deriv_order,
                        points):
    "Tabulate psi table for different integral types."
    # Handle case when list of points is empty
    if points is None:
        return _tabulate_empty_psi_table(tdim, deriv_order, fiat_element)

    # Otherwise, call FIAT to tabulate
    entity_dim = domain_to_entity_dim(integral_type, tdim)
    num_entities = num_cell_entities[cellname][entity_dim]

    mapped_points = { entity: _map_entity_points(cellname, tdim, points, entity_dim, entity)
                      for entity in range(num_entities) }

    psi_table = { entity: fiat_element.tabulate(deriv_order, entity_points)
                  for entity, entity_points in mapped_points.items() }

    return psi_table


# FIXME: Argument list here is from get_ffc_table_values, adjust to what's needed
def compute_table_values(tables, entitytype, num_points, ufl_element, flat_component, derivative_counts, epsilon):
    
    # FIXME: Instead of extracting from psi_tables, compute table here

    num_derivatives = FIXME(derivative_counts)
    points = FIXME
    integral_type = FIXME
    cellname = FIXME
    tdim = FIXME

    fiat_element = create_element(ufl_element)

    # Tabulate table of basis functions and derivatives in points
    element_table = _tabulate_psi_table(
        integral_type,
        cellname,
        tdim,
        fiat_element,
        num_derivatives,
        points)

    return element_table
'''



def missing_args(ufl_element, integral_type, cell, points):
    fiat_element = create_element(ufl_element)

    cellname = cell.cellname()
    tdim = cell.topological_dimension()

    entity_dim = domain_to_entity_dim(integral_type, tdim)
    num_entities = num_cell_entities[cellname][entity_dim]

    for entity in range(num_entities):
        entity_points = _map_entity_points(cellname, tdim, points, entity_dim, entity)
        table[entity] = fiat_element.tabulate(deriv_order, entity_points)


#tbl = compute_table(num_points, element, avg, entity, derivative_counts)
def get_ffc_table_values(psi_tables, num_points, element, avg,
                         entitytype, derivative_counts,
                         flat_component, epsilon):
    """Extract values from ffc element table.

    Returns a 3D numpy array with axes
    (entity number, quadrature point number, dof number)
    """
    # Get quadrule/element subtable
    subtable = psi_tables[num_points][element][avg]

    # Figure out shape of final array by inspecting psi_tables
    num_entities = len(subtable)
    entity = 0  # Just look at the first entity
    num_dofs = len(subtable[entity][derivative_counts])

    # Make 3D array for final result
    shape = (num_entities, num_points, num_dofs)
    res = numpy.zeros(shape)

    # Loop over entities and fill table blockwise (each block = points x dofs)
    sh = element.value_shape()
    for entity in range(num_entities):
        # Access subtable
        tbl = subtable[entity][derivative_counts]

        # Extract array for right component and order axes as (points, dofs)
        if sh == ():
            arr = numpy.transpose(tbl)
        elif len(sh) == 2 and element.num_sub_elements() == 0:
            # 2-tensor-valued elements, not a tensor product
            # mapping flat_component back to tensor component
            (_, f2t) = build_component_numbering(sh, element.symmetry())
            t_comp = f2t[flat_component]
            arr = numpy.transpose(tbl[:, t_comp[0], t_comp[1], :])
        else:
            arr = numpy.transpose(tbl[:, flat_component,:])

        # Assign block of values for this entity
        res[entity,:,:] = arr

    # Clamp almost-zeros to zero
    res[numpy.where(numpy.abs(res) < epsilon)] = 0.0
    return res


def generate_psi_table_name(element_counter,
    flat_component, derivative_counts, averaged, entitytype, num_points):
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
