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
# along with UFLACS. If not, see <http://www.gnu.org/licenses/>.

"""Tools for precomputed tables of terminal values."""

import numpy

from ufl.utils.sequences import product
from ufl.utils.derivativetuples import derivative_listing_to_counts
from ufl.permutation import build_component_numbering
from ufl.classes import FormArgument, GeometricQuantity, SpatialCoordinate, Jacobian
from ufl.algorithms.analysis import unique_tuple

from ffc.log import error

from ffc.uflacs.elementtables.table_utils import generate_psi_table_name, get_ffc_table_values
from ffc.uflacs.elementtables.table_utils import clamp_table_small_integers, strip_table_zeros, build_unique_tables

from ffc.uflacs.backends.ffc.common import ufc_restriction_offset


def get_modified_terminal_element(mt):
    gd = mt.global_derivatives
    ld = mt.local_derivatives

    # Extract element from FormArguments and relevant GeometricQuantities
    if isinstance(mt.terminal, FormArgument):
        if gd and mt.reference_value:
            error("Global derivatives of reference values not defined.")
        elif ld and not mt.reference_value:
            error("Local derivatives of global values not defined.")
        element = mt.terminal.ufl_element()
        fc = mt.flat_component
    elif isinstance(mt.terminal, SpatialCoordinate):
        if mt.reference_value:
            error("Not expecting reference value of x.")
        if gd:
            error("Not expecting global derivatives of x.")
        element = mt.terminal.ufl_domain().ufl_coordinate_element()
        if not ld:
            fc = mt.flat_component
        else:
            # Actually the Jacobian expressed as reference_grad(x)
            fc = mt.flat_component  # x-component
            assert len(mt.component) == 1
            assert mt.component[0] == mt.flat_component
    elif isinstance(mt.terminal, Jacobian):
        if mt.reference_value:
            error("Not expecting reference value of J.")
        if gd:
            error("Not expecting global derivatives of J.")
        element = mt.terminal.ufl_domain().ufl_coordinate_element()
        # Translate component J[i,d] to x element context rgrad(x[i])[d]
        assert len(mt.component) == 2
        fc, d = mt.component  # x-component, derivative
        ld = tuple(sorted((d,) + ld))
    else:
        return None

    assert not (mt.averaged and (ld or gd))

    # Change derivatives format for table lookup
    #gdim = mt.terminal.ufl_domain().geometric_dimension()
    #global_derivatives = derivative_listing_to_counts(gd, gdim)

    # Change derivatives format for table lookup
    tdim = mt.terminal.ufl_domain().topological_dimension()
    local_derivatives = derivative_listing_to_counts(ld, tdim)
    
    return element, mt.averaged, local_derivatives, fc


def build_element_tables(num_points, quadrature_rules,
                         cell, integral_type, entitytype,
                         modified_terminals, epsilon):
    """Build the element tables needed for a list of modified terminals.

    Input:
      entitytype - str
      modified_terminals - ordered sequence of unique modified terminals
      FIXME: Document

    Output:
      tables - dict(name: table)
      mt_table_names - dict(ModifiedTerminal: name)

    """
    mt_table_names = {}
    tables = {}
    table_origins = {}

    # Add to element tables
    analysis = {}
    for mt in modified_terminals:
        # FIXME: Use a namedtuple for res
        res = get_modified_terminal_element(mt)
        if res:
            analysis[mt] = res

    # Build element numbering using topological
    # ordering so subelements get priority
    from ffc.analysis import extract_sub_elements, sort_elements, _compute_element_numbers
    all_elements = [res[0] for res in analysis.values()]
    unique_elements = sort_elements(extract_sub_elements(all_elements))
    element_numbers = _compute_element_numbers(unique_elements)

    def add_table(res):
        element, avg, local_derivatives, flat_component = res

        # Build name for this particular table
        element_number = element_numbers[element]
        name = generate_psi_table_name(
            num_points, element_number, avg,
            entitytype, local_derivatives, flat_component)

        # Extract the values of the table from ffc table format
        if name not in tables:
            tables[name] = get_ffc_table_values(
                quadrature_rules[num_points][0],
                cell, integral_type,
                num_points, element, avg,
                entitytype, local_derivatives, flat_component,
                epsilon)

            # Track table origin for custom integrals:
            table_origins[name] = res
        return name

    for mt in modified_terminals:
        res = analysis.get(mt)
        if not res:
            continue
        element, avg, local_derivatives, flat_component = res

        # Generate tables for each subelement in topological ordering,
        # using same avg and local_derivatives, for each component.
        # We want the first table to be the innermost subelement so that's
        # the one the optimized tables get the name from and so that's
        # the one the table origins point to for custom integrals.
        # This results in some superfluous tables but those will be
        # removed before code generation and it's not believed to be
        # a bottleneck.
        for subelement in sort_elements(extract_sub_elements([element])):
            for fc in range(product(subelement.reference_value_shape())):
                subres = (subelement, avg, local_derivatives, fc)
                name_ignored = add_table(subres)

        # Generate table and store table name with modified terminal
        name = add_table(res)
        mt_table_names[mt] = name

    return tables, mt_table_names, table_origins


def optimize_element_tables(tables, mt_table_names, table_origins, epsilon):
    """Optimize tables and make unique set.

    Steps taken:

      - clamp values that are very close to -1, 0, +1 to those values
      - remove dofs from beginning and end of tables where values are all zero
      - for each modified terminal, provide the dof range that a given table corresponds to

    Terminology:
      name - str, name used in input arguments here
      mt - modified terminal
      table - numpy array of float values
      stripped_table - numpy array of float values with zeroes
                       removed from each end of dofrange

    Input:
      tables - { name: table }
      mt_table_names - { mt: name }

    Output:
      unique_tables - { unique_name: stripped_table }
      mt_table_ranges - { mt: (unique_name, begin, end) }
    """
    # Find and sort all unique table names mentioned in mt_table_names
    used_names = set(mt_table_names.values())
    assert None not in used_names
    #used_names.remove(None)
    used_names = sorted(used_names)

    # Drop unused tables (if any at this point)
    tables = { name: tables[name] for name in tables if name in used_names }

    # Clamp almost -1.0, 0.0, and +1.0 values first
    # (i.e. 0.999999 -> 1.0 if within epsilon distance)
    for name in used_names:
        tables[name] = clamp_table_small_integers(tables[name], epsilon)

    # Strip contiguous zero blocks at the ends of all tables
    table_ranges = {}
    for name in used_names:
        begin, end, stripped_table = strip_table_zeros(tables[name], epsilon)
        tables[name] = stripped_table
        table_ranges[name] = (begin, end)

    # Build unique table mapping
    unique_tables_list, name_to_unique_index = build_unique_tables(tables, epsilon)

    # Build mapping of constructed table names to unique names.
    # Picking first constructed name preserves some information
    # about the table origins although some names may be dropped.
    unique_names = {}
    for name in used_names:
        ui = name_to_unique_index[name]
        if ui not in unique_names:
            unique_names[ui] = name

    # Build mapping from unique table name to the table itself
    unique_tables = {}
    for ui in range(len(unique_tables_list)):
        unique_tables[unique_names[ui]] = unique_tables_list[ui]

    unique_table_origins = {}
    for ui in range(len(unique_tables_list)):
        uname = unique_names[ui]
        # Track table origins for runtime recomputation in custom integrals:
        dofrange = table_ranges[uname]
        # FIXME: Make sure the "smallest" element is chosen
        (element, avg, derivative_counts, fc) = table_origins[name]
        unique_table_origins[uname] = (element, avg, derivative_counts, fc, dofrange)

    # Build mapping from modified terminal to compacted table and dof range
    # { mt: (unique name, table dof range begin, table dof range end) }
    mt_table_ranges = {}
    for mt, name in mt_table_names.items():
        assert name is not None
        b, e = table_ranges[name]
        ui = name_to_unique_index[name]
        unique_name = unique_names[ui]
        mt_table_ranges[mt] = (unique_name, b, e)

    return unique_tables, mt_table_ranges, unique_table_origins


def offset_restricted_table_ranges(mt_table_ranges, mt_table_names,
                                   tables, modified_terminals):
    # Modify dof ranges for restricted form arguments
    # (geometry gets padded variable names instead)
    for mt in modified_terminals:
        if mt.restriction and isinstance(mt.terminal, FormArgument):
            # offset = 0 or number of dofs before table optimization
            num_original_dofs = int(tables[mt_table_names[mt]].shape[-1])
            offset = ufc_restriction_offset(mt.restriction, num_original_dofs)
            (unique_name, b, e) = mt_table_ranges[mt]
            mt_table_ranges[mt] = (unique_name, b + offset, e + offset)
    return mt_table_ranges


def analyse_table_types(unique_tables, mt_table_ranges, epsilon):
    table_types = {}
    for unique_name, table in unique_tables.items():
        #num_entities, num_points, num_dofs = table.shape
        num_points = table.shape[1]
        if product(table.shape) == 0 or numpy.allclose(table, numpy.zeros(table.shape)):  #, atol=epsilon):
            # All values are 0.0
            tabletype = "zeros"
            # All table ranges referring to this table should be empty
            assert all(data[1] == data[2]
                       for mt, data in mt_table_ranges.items()
                       if data is not None and data[0] == unique_name)
        elif numpy.allclose(table, numpy.ones(table.shape)):  #, atol=epsilon
            # All values are 1.0
            tabletype = "ones"
        elif all(numpy.allclose(table[:, 0, :], table[:, i, :])  #, atol=epsilon
                 for i in range(1, num_points)):
            # Piecewise constant over points (separately on each entity)
            tabletype = "piecewise"
        else:
            # Varying over points
            tabletype = "varying"
            # No table ranges referring to this table should be averaged
            assert all(not mt.averaged
                       for mt, data in mt_table_ranges.items()
                       if data is not None and data[0] == unique_name)
            
        table_types[unique_name] = tabletype
    return table_types


def build_optimized_tables(num_points, quadrature_rules,
                           cell, integral_type, entitytype,
                           modified_terminals, parameters):
    # Get tolerance for checking table values against 0.0 or 1.0
    from ffc.uflacs.language.format_value import get_float_threshold
    epsilon = get_float_threshold()
    # FIXME: Should be epsilon from ffc parameters
    #epsilon = parameters["epsilon"]

    # Build tables needed by all modified terminals
    tables, mt_table_names, table_origins = \
        build_element_tables(num_points, quadrature_rules,
            cell, integral_type, entitytype,
            modified_terminals, epsilon)

    # Optimize tables and get table name and dofrange for each modified terminal
    unique_tables, mt_table_ranges, table_origins = \
        optimize_element_tables(tables, mt_table_names, table_origins, epsilon)

    # Analyze tables for properties useful for optimization
    table_types = analyse_table_types(unique_tables, mt_table_ranges, epsilon)

    # Add offsets to dof ranges for restricted terminals
    mt_table_ranges = offset_restricted_table_ranges(
        mt_table_ranges, mt_table_names, tables, modified_terminals)

    # Delete unused tables and compress piecewise constant tables
    used_names = set(tabledata[0] for tabledata in mt_table_ranges.values())
    unused_names = set(unique_tables.keys()) - used_names
    for uname in unused_names:
        del table_types[uname]
        del unique_tables[uname]
    for uname, tabletype in table_types.items():
        if tabletype == "piecewise":
            # Reduce table to dimension 1 along num_points axis in generated code
            # FIXME: Make sure it's never indexed with iq!
            unique_tables[uname] = unique_tables[uname][:,0:1,:]
        elif tabletype == "zeros":
            del unique_tables[uname]
        elif tabletype == "ones":
            del unique_tables[uname]

    return unique_tables, mt_table_ranges, table_types
