# -*- coding: utf-8 -*-
# Copyright (C) 2011-2017 Martin Sandve Alnæs
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

from collections import namedtuple
import numpy

from ufl.cell import num_cell_entities
from ufl.utils.sequences import product
from ufl.utils.derivativetuples import derivative_listing_to_counts
from ufl.permutation import build_component_numbering
from ufl.classes import FormArgument, GeometricQuantity, SpatialCoordinate, Jacobian
from ufl.algorithms.analysis import unique_tuple
from ufl.measure import custom_integral_types

from ffc.log import error
from ffc.fiatinterface import create_element
from ffc.representationutils import integral_type_to_entity_dim, map_integral_points
from ffc.representationutils import create_quadrature_points_and_weights
from ffc.uflacs.backends.ffc.common import ufc_restriction_offset


# Using same defaults as numpy.allclose
default_rtol = 1e-5
default_atol = 1e-8

table_origin_t = namedtuple("table_origin",
    ["element", "avg", "derivatives", "flat_component", "dofrange", "dofmap"])


piecewise_ttypes = ("piecewise", "fixed", "ones", "zeros")


uniform_ttypes = ("uniform", "fixed", "ones", "zeros")


valid_ttypes = set(("quadrature",)) | set(piecewise_ttypes) | set(uniform_ttypes)


unique_table_reference_t = namedtuple("unique_table_reference",
    ["name", "values",
     "dofrange", "dofmap", "original_dim",
     "ttype", "is_piecewise", "is_uniform"])


def equal_tables(a, b, rtol=default_rtol, atol=default_atol):
    a = numpy.asarray(a)
    b = numpy.asarray(b)
    if a.shape != b.shape:
        return False
    else:
        return numpy.allclose(a, b, rtol=rtol, atol=atol)


def clamp_table_small_numbers(table, rtol=default_rtol, atol=default_atol, numbers=(-1.0, -0.5, 0.0, 0.5, 1.0)):
    "Clamp almost 0,1,-1 values to integers. Returns new table."
    # Get shape of table and number of columns, defined as the last axis
    table = numpy.asarray(table)
    for n in numbers:
        table[numpy.where(numpy.isclose(table, n, rtol=rtol, atol=atol))] = n
    return table


def strip_table_zeros(table, compress_zeros, rtol=default_rtol, atol=default_atol):
    "Strip zero columns from table. Returns column range (begin, end) and the new compact table."
    # Get shape of table and number of columns, defined as the last axis
    table = numpy.asarray(table)
    sh = table.shape

    # Find nonzero columns
    z = numpy.zeros(sh[:-1])  # Correctly shaped zero table
    dofmap = tuple(i for i in range(sh[-1])
                   if not numpy.allclose(z, table[..., i], rtol=rtol, atol=atol))
    if dofmap:
        # Find first nonzero column
        begin = dofmap[0]
        # Find (one beyond) last nonzero column
        end = dofmap[-1] + 1
    else:
        begin = 0
        end = 0

    # If compression is not wanted, pretend whole range is nonzero
    if not compress_zeros:
        dofmap = tuple(range(begin, end))

    # Make subtable by dropping zero columns
    stripped_table = table[..., dofmap]
    dofrange = (begin, end)
    return dofrange, dofmap, stripped_table


def build_unique_tables(tables, rtol=default_rtol, atol=default_atol):
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
            if equal_tables(u, t, rtol=rtol, atol=atol):
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
                         flat_component):
    """Extract values from ffc element table.

    Returns a 3D numpy array with axes
    (entity number, quadrature point number, dof number)
    """
    deriv_order = sum(derivative_counts)

    if integral_type in custom_integral_types:
        # Use quadrature points on cell for analysis in custom integral types
        integral_type = "cell"
        assert not avg

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
                         modified_terminals, rtol=default_rtol, atol=default_atol):
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
                entitytype, local_derivatives, flat_component)

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


def optimize_element_tables(tables, table_origins, compress_zeros, rtol=default_rtol, atol=default_atol):
    """Optimize tables and make unique set.

    Steps taken:

      - clamp values that are very close to -1, 0, +1 to those values
      - remove dofs from beginning and end of tables where values are all zero
      - for each modified terminal, provide the dof range that a given table corresponds to

    Terminology:
      name - str, name used in input arguments here
      table - numpy array of float values
      stripped_table - numpy array of float values with zeroes
                       removed from each end of dofrange

    Input:
      tables - { name: table }
      table_origins - FIXME

    Output:
      unique_tables - { unique_name: stripped_table }
      unique_table_origins - FIXME
    """
    used_names = sorted(tables)
    compressed_tables = {}
    table_ranges = {}
    table_dofmaps = {}
    table_original_num_dofs = {}

    for name in used_names:
        tbl = tables[name]

        # Clamp to selected small numbers if close,
        # (-1.0, -0.5, 0.0, 0.5 and +1.0)
        # (i.e. 0.999999 -> 1.0 if within rtol/atol distance)
        tbl = clamp_table_small_numbers(tbl, rtol=rtol, atol=atol)

        # Store original dof dimension before compressing
        num_dofs = tbl.shape[2]

        # Strip contiguous zero blocks at the ends of all tables
        dofrange, dofmap, tbl = strip_table_zeros(tbl, compress_zeros, rtol=rtol, atol=atol)

        compressed_tables[name] = tbl
        table_ranges[name] = dofrange
        table_dofmaps[name] = dofmap
        table_original_num_dofs[name] = num_dofs

    # Build unique table mapping
    unique_tables_list, name_to_unique_index = build_unique_tables(
        compressed_tables, rtol=rtol, atol=atol)

    # Build mapping of constructed table names to unique names.
    # Picking first constructed name preserves some information
    # about the table origins although some names may be dropped.
    unique_names = {}
    for name in used_names:
        ui = name_to_unique_index[name]
        if ui not in unique_names:
            unique_names[ui] = name
    table_unames = { name: unique_names[name_to_unique_index[name]]
                     for name in name_to_unique_index }

    # Build mapping from unique table name to the table itself
    unique_tables = {}
    for ui, tbl in enumerate(unique_tables_list):
        uname = unique_names[ui]
        unique_tables[uname] = tbl

    unique_table_origins = {}
    #for ui in range(len(unique_tables_list)):
    for ui in []:  # FIXME
        uname = unique_names[ui]

        # Track table origins for runtime recomputation in custom integrals:
        name = "name of 'smallest' element we can use to compute this table"
        dofrange = "FIXME"  # table_ranges[name]
        dofmap = "FIXME"  # table_dofmaps[name]

        # FIXME: Make sure the "smallest" element is chosen
        (element, avg, derivative_counts, fc) = table_origins[name]
        unique_table_origins[uname] = table_origin_t(element, avg, derivative_counts, fc, dofrange, dofmap)

    return unique_tables, unique_table_origins, table_unames, table_ranges, table_dofmaps, table_original_num_dofs


def is_zeros_table(table, rtol=default_rtol, atol=default_atol):
    return (product(table.shape) == 0
            or numpy.allclose(table, numpy.zeros(table.shape), rtol=rtol, atol=atol))


def is_ones_table(table, rtol=default_rtol, atol=default_atol):
    return numpy.allclose(table, numpy.ones(table.shape), rtol=rtol, atol=atol)


def is_quadrature_table(table, rtol=default_rtol, atol=default_atol):
    num_entities, num_points, num_dofs = table.shape
    Id = numpy.eye(num_points)
    return (num_points == num_dofs
            and all(numpy.allclose(table[i, :, :], Id, rtol=rtol, atol=atol)
                    for i in range(num_entities)))


def is_piecewise_table(table, rtol=default_rtol, atol=default_atol):
    return all(numpy.allclose(table[:, 0, :], table[:, i, :], rtol=rtol, atol=atol)
               for i in range(1, table.shape[1]))


def is_uniform_table(table, rtol=default_rtol, atol=default_atol):
    return all(numpy.allclose(table[0, :, :], table[i, :, :], rtol=rtol, atol=atol)
               for i in range(1, table.shape[0]))


def analyse_table_type(table, rtol=default_rtol, atol=default_atol):
    num_entities, num_points, num_dofs = table.shape
    if is_zeros_table(table, rtol=rtol, atol=atol):
        # Table is empty or all values are 0.0
        ttype = "zeros"
    elif is_ones_table(table, rtol=rtol, atol=atol):
        # All values are 1.0
        ttype = "ones"
    elif is_quadrature_table(table, rtol=rtol, atol=atol):
        # Identity matrix mapping points to dofs (separately on each entity)
        ttype = "quadrature"
    else:
        # Equal for all points on a given entity
        piecewise = is_piecewise_table(table, rtol=rtol, atol=atol)

        # Equal for all entities
        uniform = is_uniform_table(table, rtol=rtol, atol=atol)

        if piecewise and uniform:
            # Constant for all points and all entities
            ttype = "fixed"
        elif piecewise:
            # Constant for all points on each entity separately
            ttype = "piecewise"
        elif uniform:
            # Equal on all entities
            ttype = "uniform"
        else:
            # Varying over points and entities
            ttype = "varying"
    return ttype


def analyse_table_types(unique_tables, rtol=default_rtol, atol=default_atol):
    return { uname: analyse_table_type(table, rtol=rtol, atol=atol)
             for uname, table in unique_tables.items() }


def build_optimized_tables(num_points, quadrature_rules,
                           cell, integral_type, entitytype,
                           modified_terminals, existing_tables,
                           compress_zeros, rtol=default_rtol, atol=default_atol):
    # Build tables needed by all modified terminals
    tables, mt_table_names, table_origins = \
        build_element_tables(num_points, quadrature_rules,
            cell, integral_type, entitytype,
            modified_terminals, rtol=rtol, atol=atol)

    # Optimize tables and get table name and dofrange for each modified terminal
    unique_tables, unique_table_origins, table_unames, table_ranges, table_dofmaps, table_original_num_dofs = \
        optimize_element_tables(tables, table_origins, compress_zeros, rtol=rtol, atol=atol)

    # Get num_dofs for all tables before they can be deleted later
    unique_table_num_dofs = { uname: tbl.shape[2] for uname, tbl in unique_tables.items() }

    # Analyze tables for properties useful for optimization
    unique_table_ttypes = analyse_table_types(unique_tables, rtol=rtol, atol=atol)

    # Compress tables that are constant along num_entities or num_points
    for uname, tabletype in unique_table_ttypes.items():
        if tabletype in piecewise_ttypes:
            # Reduce table to dimension 1 along num_points axis in generated code
            unique_tables[uname] = unique_tables[uname][:,0:1,:]
        if tabletype in uniform_ttypes:
            # Reduce table to dimension 1 along num_entities axis in generated code
            unique_tables[uname] = unique_tables[uname][0:1,:,:]

    # Delete tables not referenced by modified terminals
    used_unames = set(table_unames[name] for name in mt_table_names.values())
    unused_unames = set(unique_tables.keys()) - used_unames
    for uname in unused_unames:
        del unique_table_ttypes[uname]
        del unique_tables[uname]

    # Change tables to point to existing optimized tables
    # (i.e. tables from other contexts that have been compressed to look the same)
    name_map = {}
    existing_names = sorted(existing_tables)
    for uname in sorted(unique_tables):
        utbl = unique_tables[uname]
        for i, ename in enumerate(existing_names):
            etbl = existing_tables[ename]
            if equal_tables(utbl, etbl, rtol=rtol, atol=atol):
                # Setup table name mapping
                name_map[uname] = ename
                # Don't visit this table again (just to avoid the processing)
                existing_names.pop(i)
                break

    # Replace unique table names
    for uname, ename in name_map.items():
        unique_tables[ename] = existing_tables[ename]
        del unique_tables[uname]
        unique_table_ttypes[ename] = unique_table_ttypes[uname]
        del unique_table_ttypes[uname]

    # Build mapping from modified terminal to unique table with metadata
    # { mt: (unique name,
    #        (table dof range begin, table dof range end),
    #        [top parent element dof index for each local index],
    #        ttype, original_element_dim) }
    mt_unique_table_reference = {}
    for mt, name in list(mt_table_names.items()):
        # Get metadata for the original table (name is not the unique name!)
        dofrange = table_ranges[name]
        dofmap = table_dofmaps[name]
        original_dim = table_original_num_dofs[name]

        # Map name -> uname
        uname = table_unames[name]

        # Map uname -> ename
        ename = name_map.get(uname, uname)

        # Some more metadata stored under the ename
        ttype = unique_table_ttypes[ename]

        # Add offset to dofmap and dofrange for restricted terminals
        if mt.restriction and isinstance(mt.terminal, FormArgument):
            # offset = 0 or number of dofs before table optimization
            offset = ufc_restriction_offset(mt.restriction, original_dim)
            (b, e) = dofrange
            dofrange = (b + offset, e + offset)
            dofmap = tuple(i + offset for i in dofmap)

        # Store reference to unique table for this mt
        mt_unique_table_reference[mt] = unique_table_reference_t(
            ename, unique_tables[ename],
            dofrange, dofmap, original_dim,
            ttype, ttype in piecewise_ttypes, ttype in uniform_ttypes)

    return unique_tables, unique_table_ttypes, unique_table_num_dofs, mt_unique_table_reference
