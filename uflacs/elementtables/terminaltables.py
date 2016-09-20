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
# along with UFLACS. If not, see <http://www.gnu.org/licenses/>.

"""Tools for precomputed tables of terminal values."""

import numpy

from ufl.utils.sequences import product
from ufl.utils.derivativetuples import derivative_listing_to_counts
from ufl.classes import FormArgument, GeometricQuantity, SpatialCoordinate, Jacobian
from ufl.algorithms.analysis import unique_tuple

from ffc.log import error

from uflacs.elementtables.table_utils import generate_psi_table_name, get_ffc_table_values
from uflacs.elementtables.table_utils import clamp_table_small_integers, strip_table_zeros, build_unique_tables

from uflacs.backends.ffc.common import ufc_restriction_offset


def build_element_tables(psi_tables, num_points, entitytype, modified_terminals, epsilon):
    """Build the element tables needed for a list of modified terminals.

    Input:
      psi_tables - tables from ffc
      entitytype - str
      modified_terminals - ordered sequence of unique modified terminals

    Output:
      tables - dict(name: table)
      mt_table_names - dict(ModifiedTerminal: name)

    """
    element_counter_map = {}
    mt_table_names = {}
    tables = {}
    # Add to element tables
    for mt in modified_terminals:
        
        t = mt.terminal
        rv = mt.reference_value
        gd = mt.global_derivatives
        ld = mt.local_derivatives
        gc = mt.component
        fc = mt.flat_component

        # Extract element from FormArguments and relevant GeometricQuantities
        if isinstance(t, FormArgument):
            if gd and rv:
                error("Global derivatives of reference values not defined.")
            elif ld and not rv:
                error("Local derivatives of global values not defined.")
            element = t.ufl_element()
        elif isinstance(t, SpatialCoordinate):
            if rv:
                error("Not expecting reference value of x.")
            if gd:
                error("Not expecting global derivatives of x.")
            element = t.ufl_domain().ufl_coordinate_element()
            if ld:
                # Actually the Jacobian, translate component gc to x element context
                fc, ld = gc
                ld = (ld,)
        elif isinstance(t, Jacobian):
            if rv:
                error("Not expecting reference value of J.")
            if gd:
                error("Not expecting global derivatives of J.")
            element = t.ufl_domain().ufl_coordinate_element()
            fc = gc[0]
            ld = tuple(sorted((gc[1],) + ld))
        else:
            continue

        # Count elements as we go
        element_counter = element_counter_map.get(element)
        if element_counter is None:
            element_counter = len(element_counter_map)
            element_counter_map[element] = element_counter

        # Change derivatives format for table lookup
        #if gd:
        #    gdim = t.ufl_domain().geometric_dimension()
        #    global_derivatives = tuple(derivative_listing_to_counts(gd, gdim))
        #else:
        #    global_derivatives = None

        # Change derivatives format for table lookup
        if ld:
            tdim = t.ufl_domain().topological_dimension()
            local_derivatives = tuple(derivative_listing_to_counts(ld, tdim))
        else:
            local_derivatives = None

        # Build name for this particular table
        name = generate_psi_table_name(element_counter, fc, local_derivatives,
                                       mt.averaged, entitytype, num_points)

        # Extract the values of the table from ffc table format
        table = tables.get(name)
        if table is None:
            table = get_ffc_table_values(psi_tables, entitytype, num_points,
                                         element, fc, local_derivatives, epsilon)
            tables[name] = table

        # Store table name with modified terminal
        mt_table_names[mt] = name

    return tables, mt_table_names


def optimize_element_tables(tables, mt_table_names, epsilon):
    """Optimize tables and make unique set.

    Steps taken:
    - clamp values that are very close to -1, 0, +1 to those values
    - remove dofs from beginning and end of tables where values are all zero
    - for each modified terminal, provide the dof range that a given table corresponds to

    Input:
      tables - a mapping from name to table values
      mt_table_names - a mapping from modified terminal to table name

    Output:
      unique_tables - a mapping from name to table values with stripped zero columns
      mt_table_ranges - a mapping from modified terminal to (name, begin, end)
    """
    # Names here are a bit long and slightly messy...

    # Drop tables not mentioned in mt_table_names
    used_names = set(mt_table_names.values())
    assert None not in mt_table_names
    #used_names.remove(None)
    tables = { name: tables[name]
               for name in tables
               if name in used_names }

    # Clamp almost -1.0, 0.0, and +1.0 values first
    # (i.e. 0.999999 -> 1.0 if within epsilon distance)
    tables = { name: clamp_table_small_integers(table, epsilon)
               for name, table in tables.items() }

    # Strip contiguous zero blocks at the ends of all tables
    stripped_tables = {}
    table_ranges = {}
    for name, table in tables.items():
        begin, end, stripped_table = strip_table_zeros(table, epsilon)
        # Drop empty tables
        if product(stripped_table.shape) == 0:
            end = begin
        else:
            stripped_tables[name] = stripped_table
        table_ranges[name] = (begin, end)

    # Build unique table mapping
    unique_tables_list, table_name_to_unique_index = build_unique_tables(stripped_tables, epsilon)

    # Build mapping of constructed table names to unique names.
    # Picking first constructed name preserves some information
    # about the table origins although some names may be dropped.
    unique_table_names = {}
    for name in sorted(table_name_to_unique_index):
        unique_index = table_name_to_unique_index[name]
        if unique_index in unique_table_names:
            continue
        unique_table_names[unique_index] = name

    # Build mapping from unique table name to the table itself
    unique_tables = { unique_table_names[unique_index]: unique_tables_list[unique_index]
                      for unique_index in range(len(unique_tables_list)) }

    # Build mapping from modified terminal to compacted table data:
    # mt ->
    #    (unique name, table dof range begin, table dof range end) for varying tables
    #    (None, begin, end=begin) for empty range tables
    #    TODO:  ("ones", begin, end) for tables where all values are 1.0
    #    TODO:  ("zeros", begin, end) for tables where all values are 0.0
    mt_table_ranges = {}
    for mt, name in mt_table_names.items():
        if name is not None:
            b, e = table_ranges[name]
            if e - b > 0:
                unique_index = table_name_to_unique_index[name]
                unique_name = unique_table_names[unique_index]
            else:
                unique_name = None
            mt_table_ranges[mt] = (unique_name, b, e)

    return unique_tables, mt_table_ranges


class TableProvider(object):
    def __init__(self, psi_tables, parameters):
        self.psi_tables = psi_tables

        # FIXME: Should be epsilon from ffc parameters
        from uflacs.language.format_value import get_float_threshold
        self.epsilon = get_float_threshold()

    def build_optimized_tables(self, num_points, entitytype, modified_terminals):
        psi_tables = self.psi_tables
        epsilon = self.epsilon

        # FIXME: Refactor such that
        #        terminal_table_ranges = dict(modified_terminal -> (tablename, begin, end))
        #    instead of list where modified terminal indices have meaning,
        #    to avoid having to renumber later! Requires modified terminals to be hashable, are they?
        #    Note: Include modification for restricted form arguments
        #    Note: Do not include DG0 and Real component tables (or use tablename "ones" which can be used to define piecewise constants in partitioning and also checked for in code generation)
        #    Note: Do not include empty tables (why do they occur?) (or use tablename "empty" and begin=end)
        #    Or split terminal_table_ranges into (varying_table_ranges, constant_tables, empty_tables)?
        #    If the code generation

        # Build tables needed by all modified terminals
        # (currently build here means extract from ffc psi_tables)
        tables, mt_table_names = \
            build_element_tables(psi_tables, num_points, entitytype, modified_terminals, epsilon)

        # Optimize tables and get table name and dofrange for each modified terminal
        unique_tables, mt_table_ranges = \
            optimize_element_tables(tables, mt_table_names, epsilon)

        # Modify dof ranges for restricted form arguments
        # (geometry gets padded variable names instead)
        for mt in modified_terminals:
            if mt.restriction and isinstance(mt.terminal, FormArgument):
                # offset = 0 or number of dofs before table optimization
                num_original_dofs = int(tables[mt_table_names[mt]].shape[-1])
                offset = ufc_restriction_offset(mt.restriction, num_original_dofs)
                (unique_name, b, e) = mt_table_ranges[mt]
                mt_table_ranges[mt] = (unique_name, b + offset, e + offset)

        return unique_tables, mt_table_ranges
