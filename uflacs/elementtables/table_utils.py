# -*- coding: utf-8 -*-
# Copyright (C) 2011-2015 Martin Sandve Alnæs
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

from __future__ import print_function # used in some debugging

from six import itervalues, iterkeys
from six import advance_iterator as next
from six.moves import map
from six.moves import xrange as range
import numpy as np

default_tolerance = 1e-14


def equal_tables(a, b, eps=default_tolerance):
    "Compare tables to be equal within a tolerance."
    a = np.asarray(a)
    b = np.asarray(b)
    if a.shape != b.shape:
        return False
    if len(a.shape) > 1:
        return all(equal_tables(a[i], b[i], eps) for i in range(a.shape[0]))

    def scalars_equal(x, y, eps):
        return abs(x-y) < eps
    return all(scalars_equal(a[i], b[i], eps) for i in range(a.shape[0]))


def strip_table_zeros(table, eps=default_tolerance):
    "Strip zero columns from table. Returns column range (begin,end) and the new compact table."
    # Get shape of table and number of columns, defined as the last axis
    table = np.asarray(table)
    sh = table.shape
    nc = sh[-1]

    # Find first nonzero column
    begin = nc
    for i in range(nc):
        if np.linalg.norm(table[..., i]) > eps:
            begin = i
            break

    # Find (one beyond) last nonzero column
    end = begin
    for i in range(nc-1, begin-1, -1):
        if np.linalg.norm(table[..., i]) > eps:
            end = i+1
            break

    # Make subtable by stripping first and last columns
    return begin, end, table[..., begin:end]


def build_unique_tables(tables):
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
            if equal_tables(u, t):
                found = i
                break
        if found == -1:
            i = len(unique)
            unique.append(t)
        mapping[k] = i
    return unique, mapping


def build_unique_tables2(tables):
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
            if equal_tables(u, t):
                found = i
                break
        if found == -1:
            i = len(unique)
            unique.append(t)
        mapping[k] = i

    return unique, mapping


def get_ffc_table_values(tables, entitytype, num_points, element, flat_component, derivative_counts):
    """Extract values from ffc element table.

    Returns a 3D numpy array with axes
    (entity number, quadrature point number, dof number)
    """
    # Get quadrule/element subtable
    element_table = tables[num_points][element]

    # Temporary fix for new table structure TODO: Handle avg properly
    if len(element_table) != 1:
        print()
        print(element_table)
    assert len(element_table) == 1
    element_table = element_table[None]

    # FFC property:
    # element_counter = element_map[num_points][element]

    # Figure out shape of final array by inspecting tables
    num_entities = len(element_table)
    tmp = next(itervalues(element_table)) # Pick subtable for arbitrary chosen cell entity
    if derivative_counts is None: # Workaround for None vs (0,)*tdim
        dc = next(iterkeys(tmp))
        derivative_counts = (0,)*len(dc)
    num_dofs = len(tmp[derivative_counts])

    # Make 3D array for final result
    shape = (num_entities, num_points, num_dofs)
    res = np.zeros(shape)

    # Loop over entities and fill table blockwise (each block = points x dofs)
    sh = element.value_shape()
    for entity in range(num_entities):
        # Access subtable
        entity_key = None if entitytype == "cell" else entity
        tbl = element_table[entity_key][derivative_counts]

        # Extract array for right component and order axes as (points, dofs)
        if sh == ():
            arr = np.transpose(tbl)
        else:
            arr = np.transpose(tbl[:, flat_component,:])

        # Assign block of values for this entity
        res[entity,:,:] = arr

    # Clamp almost-zeros to zero  # FIXME: Use parameters["epsilon"] from ffc for this tolerance
    zero_clamp_tolerance = default_tolerance
    res[np.where(np.abs(res) < zero_clamp_tolerance)] = 0.0
    return res


def generate_psi_table_name(element_counter, flat_component, derivative_counts, averaged, entitytype, num_quadrature_points):
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

    if isinstance(flat_component, int):
        name += "_C%d" % flat_component
    else:
        assert flat_component is None

    if derivative_counts and any(derivative_counts):
        name += "_D" + "".join(map(str, derivative_counts))

    if averaged == "cell":
        name += "_AC"
    elif averaged == "facet":
        name += "_AF"

    if entitytype == "cell":
        pass
    elif entitytype == "facet":
        name += "_F"
    elif entitytype == "vertex":
        name += "_V"
    else:
        error("Unknown entity type %s." % entitytype)

    if isinstance(num_quadrature_points, int):
        name += "_Q%d" % num_quadrature_points
    else:
        assert num_quadrature_points is None

    return name


def _examples(tables):
    name = generate_psi_table_name(counter, flat_component, derivative_counts, averaged, entitytype, None)
    values = get_ffc_table_values(tables, entitytype, num_points, element, flat_component, derivative_counts)

    begin, end, table = strip_table_zeros(table)
    all_zeros = table.shape[-1] == 0
    all_ones = equal_tables(table, np.ones(table.shape))
