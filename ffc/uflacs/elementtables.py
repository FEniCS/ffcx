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

from __future__ import print_function  # used in some debugging

import numpy

from ufl.cell import num_cell_entities
from ufl.utils.sequences import product
from ufl.utils.derivativetuples import derivative_listing_to_counts
from ufl.permutation import build_component_numbering
from ufl.classes import FormArgument, GeometricQuantity, SpatialCoordinate, Jacobian
from ufl.algorithms.analysis import unique_tuple

from ffc.log import error
from ffc.fiatinterface import create_element
from ffc.representationutils import integral_type_to_entity_dim, map_integral_points
from ffc.representationutils import create_quadrature_points_and_weights
from ffc.uflacs.backends.ffc.common import ufc_restriction_offset


# TODO: Use this class for tables with metadata?
class Table(object):
    """Table with metadata.

    Valid table types:
    "zeros"
    "ones"
    "quadrature"
    "piecewise"
    "uniform"
    "fixed"
    "varying"

    FIXME: Document these. For now see table computation.
    """
    def __init__(self, name, values, tabletype):
        self.name = name
        self.values = values
        self.num_entities = values.shape[0]
        self.num_points = values.shape[1]
        self.num_dofs = values.shape[2]
        self.tabletype = tabletype

        self.piecewise = tabletype in ("piecewise", "fixed")
        self.uniform = tabletype in ("uniform", "fixed")


class TableData(object):
    """

    Valid table types:
    "zeros"
    "ones"
    "quadrature"
    "piecewise"
    "uniform"
    "fixed"
    "varying"

    FIXME: Document these. For now see table computation.
    """
    def __init__(self, name, values, begin, end, tabletype):
        self.name = name
        self.begin = begin
        self.end = end
        self.values = values

        self.tabletype = tabletype
        self.is_piecewise = tabletype in ("piecewise", "fixed")
        self.is_uniform = tabletype in ("uniform", "fixed")


# TODO: Replace with numpy.allclose
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


# TODO: Is clamping 1's really safe?
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
                element, avg,
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


def is_zeros_table(table, epsilon):
    return (product(table.shape) == 0
            or numpy.allclose(table, numpy.zeros(table.shape), atol=epsilon))


def is_ones_table(table, epsilon):
    return numpy.allclose(table, numpy.ones(table.shape), atol=epsilon)


def is_quadrature_table(table, epsilon):
    num_entities, num_points, num_dofs = table.shape
    I = numpy.eye(num_points)
    return (num_points == num_dofs
            and all(numpy.allclose(table[i, :, :], I, atol=epsilon)
                    for i in range(num_entities)))


def is_piecewise_table(table, epsilon):
    return all(numpy.allclose(table[:, 0, :], table[:, i, :], atol=epsilon)
               for i in range(1, table.shape[1]))


def is_uniform_table(table, epsilon):
    return all(numpy.allclose(table[0, :, :], table[i, :, :], atol=epsilon)
               for i in range(1, table.shape[0]))


def analyse_table_types(unique_tables, epsilon):
    table_types = {}
    for unique_name, table in unique_tables.items():
        num_entities, num_points, num_dofs = table.shape
        if is_zeros_table(table, epsilon):
            # Table is empty or all values are 0.0
            tabletype = "zeros"
        elif is_ones_table(table, epsilon):
            # All values are 1.0
            tabletype = "ones"
        elif is_quadrature_table(table, epsilon):
            # Identity matrix mapping points to dofs (separately on each entity)
            tabletype = "quadrature"
        else:
            # Equal for all points on a given entity
            piecewise = is_piecewise_table(table, epsilon)

            # Equal for all entities
            uniform = is_uniform_table(table, epsilon)

            if piecewise and uniform:
                # Constant for all points and all entities
                tabletype = "fixed"
            elif piecewise:
                # Constant for all points on each entity separately
                tabletype = "piecewise"
            elif uniform:
                # Equal on all entities
                tabletype = "uniform"
            else:
                # Varying over points and entities
                tabletype = "varying"

        table_types[unique_name] = tabletype

    return table_types


def build_optimized_tables(num_points, quadrature_rules,
                           cell, integral_type, entitytype,
                           modified_terminals, existing_tables,
                           parameters):
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
    table_types = analyse_table_types(unique_tables, epsilon)

    # Get num_dofs for all tables before they can be deleted later
    table_num_dofs = { name: tbl.shape[2] for name, tbl in unique_tables.items() }

    # Consistency checking
    for unique_name, tabletype in table_types.items():
        if tabletype == "zeros":
            # All table ranges referring to this table should be empty
            assert all(data[1] == data[2]
                       for mt, data in mt_table_ranges.items()
                       if data is not None and data[0] == unique_name)
        if tabletype == "varying":
            # No table ranges referring to this table should be averaged
            assert all(not mt.averaged
                       for mt, data in mt_table_ranges.items()
                       if data is not None and data[0] == unique_name)

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
        if tabletype in ("piecewise", "fixed"):
            # Reduce table to dimension 1 along num_points axis in generated code
            unique_tables[uname] = unique_tables[uname][:,0:1,:]
        if tabletype in ("uniform", "fixed"):
            # Reduce table to dimension 1 along num_entities axis in generated code
            unique_tables[uname] = unique_tables[uname][0:1,:,:]
        if tabletype in ("zeros", "ones", "quadrature"):
            del unique_tables[uname]

    # Change tables to point to existing optimized tables
    name_map = {}
    existing_names = set(existing_tables)
    for name1 in sorted(unique_tables):
        tbl1 = unique_tables[name1]
        for name2 in existing_names:
            tbl2 = existing_tables[name2]
            if equal_tables(tbl1, tbl2, epsilon):
                # Don't visit this table again (just to avoid the processing)
                existing_names.remove(name2)
                # Replace table name
                name_map[name1] = name2
                del unique_tables[name1]
                unique_tables[name2] = tbl2
                table_types[name2] = table_types[name1]
                del table_types[name1]
                break
    for mt in list(mt_table_ranges):
        tr = mt_table_ranges[mt]
        if tr[0] in name_map:
            mt_table_ranges[mt] = (name_map[tr[0]],) + tuple(tr[1:])

    return unique_tables, mt_table_ranges, table_types, table_num_dofs
