# -*- coding: utf-8 -*-
"Utility functions for quadrature representation."

# Copyright (C) 2007-2010 Kristian B. Oelgaard
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
# First added:  2007-03-16
# Last changed: 2015-03-27
#
# Hacked by Marie E. Rognes 2013
# Modified by Anders Logg 2014
# Modified by Lizao Li 2015

# Python modules.
import numpy

# FFC modules.
from ffc.log import debug, error
from ffc.quadrature.cpp import format


def create_psi_tables(tables, eliminate_zeros, entity_type):
    "Create names and maps for tables and non-zero entries if appropriate."

    debug("\nQG-utils, psi_tables:\n" + str(tables))

    # Create element map {points:{element:number,},} and a plain
    # dictionary {name:values,}.
    element_map, flat_tables = flatten_psi_tables(tables, entity_type)
    debug("\nQG-utils, psi_tables, flat_tables:\n" + str(flat_tables))

    # Reduce tables such that we only have those tables left with
    # unique values. Create a name map for those tables that are
    # redundant.
    name_map, unique_tables = unique_psi_tables(flat_tables, eliminate_zeros)
    debug("\nQG-utils, psi_tables, unique_tables:\n" + str(unique_tables))
    debug("\nQG-utils, psi_tables, name_map:\n" + str(name_map))

    return (element_map, name_map, unique_tables)


def flatten_psi_tables(tables, entity_type):
    """Create a 'flat' dictionary of tables with unique names and a name
    map that maps number of quadrature points and element name to a
    unique element number.

    Input tables on the format for scalar and non-scalar elements respectively:
      tables[num_points][element][entity][derivs][ip][dof]
      tables[num_points][element][entity][derivs][ip][component][dof]

    Planning to change this into:
      tables[num_points][element][avg][entity][derivs][ip][dof]
      tables[num_points][element][avg][entity][derivs][ip][component][dof]

    Returns:
      element_map - { num_quad_points: {ufl_element: element_number} }.
      flat_tables - { unique_table_name: values[ip,dof] }.

    """

    generate_psi_name = format["psi name"]

    def sorted_items(mapping, **sorted_args):
        return [(k, mapping[k]) for k in sorted(mapping.keys(), **sorted_args)]

    # Initialise return values and element counter.
    flat_tables = {}
    element_map = {}
    counter = 0

    # There's a separate set of tables for each number of quadrature points
    for num_points, element_tables in sorted_items(tables):
        element_map[num_points] = {}

        # There's a set of tables for each element
        for element, avg_tables in sorted_items(element_tables,
                                                key=lambda x: str(x)):
            element_map[num_points][element] = counter

            # There's a set of tables for non-averaged and averaged
            # (averaged only occurs with num_points == 1)
            for avg, entity_tables in sorted_items(avg_tables):

                # There's a set of tables for each entity number (only
                # 1 for the cell, >1 for facets and vertices)
                for entity, derivs_tables in sorted_items(entity_tables):

                    # There's a set of tables for each derivative combination
                    for derivs, fiat_tables in sorted_items(derivs_tables):
                        # Flatten fiat_table for tensor-valued basis
                        # This is necessary for basic
                        # (non-tensor-product) tensor elements

                        if len(numpy.shape(fiat_tables)) > 3:
                            shape = fiat_tables.shape
                            value_shape = shape[1:-1]
                            fiat_tables = fiat_tables.reshape((shape[0],
                                                               numpy.product(value_shape),
                                                               shape[-1]))

                        # Transform fiat_tables to a list of tables on
                        # the form psi_table[dof][ip] for each scalar
                        # component
                        if element.value_shape():
                            # Table is on the form
                            # fiat_tables[ip][component][dof].
                            transposed_table = numpy.transpose(fiat_tables, (1, 2, 0))
                            component_tables = list(enumerate(transposed_table))
                        else:
                            # Scalar element, table is on the form
                            # fiat_tables[ip][dof].  Using () for the
                            # component because generate_psi_name
                            # expects that
                            component_tables = [((),
                                                 numpy.transpose(fiat_tables))]

                        # Iterate over the innermost tables for each
                        # scalar component
                        for component, psi_table in component_tables:

                            # Generate the table name.
                            name = generate_psi_name(counter, entity_type,
                                                     entity, component,
                                                     derivs, avg)

                            # Verify uniqueness of names
                            if name in flat_tables:
                                error("Table name is not unique, something is wrong:\n  name = %s\n  table = %s\n" % (name, flat_tables))

                            # Store table with unique name
                            flat_tables[name] = psi_table

            # Increase unique numpoints*element counter
            counter += 1

    return (element_map, flat_tables)


def unique_psi_tables(tables, eliminate_zeros):
    """Returns a name map and a dictionary of unique tables. The function
    checks if values in the tables are equal, if this is the case it
    creates a name mapping. It also create additional information
    (depending on which parameters are set) such as if the table
    contains all ones, or only zeros, and a list on non-zero columns.
    unique_tables - {name:values,}.  name_map -
    {original_name:[new_name, non-zero-columns (list), is zero (bool),
    is ones (bool)],}.

    """

    # Get unique tables (from old table utility).
    name_map, inverse_name_map = unique_tables(tables)

    debug("\ntables: " + str(tables))
    debug("\nname_map: " + str(name_map))
    debug("\ninv_name_map: " + str(inverse_name_map))

    # Set values to zero if they are lower than threshold.
    format_epsilon = format["epsilon"]
    for name in tables:
        # Get values.
        vals = tables[name]
        for r in range(numpy.shape(vals)[0]):
            for c in range(numpy.shape(vals)[1]):
                if abs(vals[r][c]) < format_epsilon:
                    vals[r][c] = 0
        tables[name] = vals

    # Extract the column numbers that are non-zero.  If optimisation
    # option is set counter for non-zero column arrays.
    i = 0
    non_zero_columns = {}
    if eliminate_zeros:
        for name in sorted(tables.keys()):

            # Get values.
            vals = tables[name]

            # Skip if values are missing
            if len(vals) == 0:
                continue

            # Use the first row as reference.
            non_zeros = list(vals[0].nonzero()[0])

            # If all columns in the first row are non zero, there's no
            # point in continuing.
            if len(non_zeros) == numpy.shape(vals)[1]:
                continue

            # If we only have one row (IP) we just need the nonzero
            # columns.
            if numpy.shape(vals)[0] == 1:
                if list(non_zeros):
                    non_zeros.sort()
                    non_zero_columns[name] = (i, non_zeros)

                    # Compress values.
                    tables[name] = vals[:, non_zeros]
                    i += 1

            # Check if the remaining rows are nonzero in the same
            # positions, else expand.
            else:
                for j in range(1, numpy.shape(vals)[0]):
                    # All rows must have the same non-zero columns for
                    # the optimization to work (at this stage).
                    new_non_zeros = list(vals[j].nonzero()[0])
                    if non_zeros != new_non_zeros:
                        non_zeros = non_zeros + [c for c in new_non_zeros if c not in non_zeros]
                        # If this results in all columns being
                        # non-zero, continue.
                        if len(non_zeros) == numpy.shape(vals)[1]:
                            continue

                # Only add nonzeros if it results in a reduction of
                # columns.
                if len(non_zeros) != numpy.shape(vals)[1]:
                    if list(non_zeros):
                        non_zeros.sort()
                        non_zero_columns[name] = (i, non_zeros)

                        # Compress values.
                        tables[name] = vals[:, non_zeros]
                        i += 1

    # Check if we have some zeros in the tables.
    names_zeros = contains_zeros(tables)

    # Get names of tables with all ones.
    names_ones = get_ones(tables)

    # Add non-zero column, zero and ones info to inverse_name_map (so
    # we only need to pass around one name_map to code generating
    # functions).
    for name in inverse_name_map:
        if inverse_name_map[name] in non_zero_columns:
            nzc = non_zero_columns[inverse_name_map[name]]
            zero = inverse_name_map[name] in names_zeros
            ones = inverse_name_map[name] in names_ones
            inverse_name_map[name] = [inverse_name_map[name], nzc, zero, ones]
        else:
            zero = inverse_name_map[name] in names_zeros
            ones = inverse_name_map[name] in names_ones
            inverse_name_map[name] = [inverse_name_map[name], (), zero, ones]

    # If we found non zero columns we might be able to reduce number
    # of tables further.
    if non_zero_columns:
        # Try reducing the tables. This is possible if some tables
        # have become identical as a consequence of compressing the
        # tables.  This happens with e.g., gradients of linear basis
        # FE0 = {-1,0,1}, nzc0 = [0,2]
        # FE1 = {-1,1,0}, nzc1 = [0,1]  -> FE0 = {-1,1}, nzc0 = [0,2], nzc1 = [0,1].

        # Call old utility function again.
        nm, inv_nm = unique_tables(tables)

        # Update name maps.
        for name in inverse_name_map:
            if inverse_name_map[name][0] in inv_nm:
                inverse_name_map[name][0] = inv_nm[inverse_name_map[name][0]]
        for name in nm:
            maps = nm[name]
            for m in maps:
                if name not in name_map:
                    name_map[name] = []
                if m in name_map:
                    name_map[name] += name_map[m] + [m]
                    del name_map[m]
                else:
                    name_map[name].append(m)

        # Get new names of tables with all ones (for vector
        # constants).
        names = get_ones(tables)

        # Because these tables now contain ones as a consequence of
        # compression we still need to consider the non-zero columns
        # when looking up values in coefficient arrays. The psi
        # entries can however we neglected and we don't need to
        # tabulate the values (if option is set).
        for name in names:
            if name in name_map:
                maps = name_map[name]
                for m in maps:
                    inverse_name_map[m][3] = True
            if name in inverse_name_map:
                    inverse_name_map[name][3] = True

    # Write protect info and return values
    for name in inverse_name_map:
        inverse_name_map[name] = tuple(inverse_name_map[name])

    # Note: inverse_name_map here is called name_map in
    # create_psi_tables and the quadraturetransformerbase class
    return (inverse_name_map, tables)


def unique_tables(tables):
    """Removes tables with redundant values and returns a name_map and a
    inverse_name_map. E.g.,

    tables = {a:[0,1,2], b:[0,2,3], c:[0,1,2], d:[0,1,2]}
    results in:
    tables = {a:[0,1,2], b:[0,2,3]}
    name_map = {a:[c,d]}
    inverse_name_map = {a:a, b:b, c:a, d:a}.

    """

    format_epsilon = format["epsilon"]

    name_map = {}
    inverse_name_map = {}
    names = sorted(tables.keys())
    mapped = []

    # Loop all tables to see if some are redundant.
    for i in range(len(names)):
        name0 = names[i]
        if name0 in mapped:
            continue
        val0 = numpy.array(tables[name0])

        for j in range(i + 1, len(names)):
            name1 = names[j]
            if name1 in mapped:
                continue
            val1 = numpy.array(tables[name1])

            # Check if dimensions match.
            if numpy.shape(val0) == numpy.shape(val1):
                # Check if values are the same.
                if len(val0) > 0 and abs(val0 - val1).max() < format_epsilon:
                    mapped.append(name1)
                    del tables[name1]
                    if name0 in name_map:
                        name_map[name0].append(name1)
                    else:
                        name_map[name0] = [name1]
                    # Create inverse name map.
                    inverse_name_map[name1] = name0

    # Add self.
    for name in tables:
        if name not in inverse_name_map:
            inverse_name_map[name] = name

    return (name_map, inverse_name_map)


def get_ones(tables):
    "Return names of tables for which all values are 1.0."
    f_epsilon = format["epsilon"]
    names = []
    for name in tables:
        vals = tables[name]
        if len(vals) > 0 and abs(vals - numpy.ones(numpy.shape(vals))).max() < f_epsilon:
            names.append(name)
    return names


def contains_zeros(tables):
    "Checks if any tables contains only zeros."
    f_epsilon = format["epsilon"]
    names = []
    for name in tables:
        vals = tables[name]
        if len(vals) > 0 and abs(vals).max() < f_epsilon:
            names.append(name)
    return names


def create_permutations(expr):

    # This is probably not used.
    if len(expr) == 0:
        return expr

    # Format keys and values to lists and tuples.
    if len(expr) == 1:
        new = {}
        for key, val in expr[0].items():
            if key == ():
                pass
            elif not isinstance(key[0], tuple):
                key = (key,)
            if not isinstance(val, list):
                val = [val]
            new[key] = val

        return new

    # Create permutations of two lists.
    # TODO: there could be a cleverer way of changing types of keys and vals.
    if len(expr) == 2:
        new = {}
        for key0, val0 in expr[0].items():
            if isinstance(key0[0], tuple):
                key0 = list(key0)
            if not isinstance(key0, list):
                key0 = [key0]
            if not isinstance(val0, list):
                val0 = [val0]
            for key1, val1 in expr[1].items():
                if key1 == ():
                    key1 = []
                elif isinstance(key1[0], tuple):
                    key1 = list(key1)
                if not isinstance(key1, list):
                    key1 = [key1]
                if not isinstance(val1, list):
                    val1 = [val1]
                if tuple(key0 + key1) in new:
                    error("This is not supposed to happen.")
                new[tuple(key0 + key1)] = val0 + val1

        return new

    # Create permutations by calling this function recursively.  This
    # is only used for rank > 2 tensors I think.
    # if len(expr) > 2:
    #     new = permutations(expr[0:2])
    #     return permutations(new + expr[2:])
