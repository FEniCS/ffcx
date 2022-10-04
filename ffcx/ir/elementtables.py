# Copyright (C) 2013-2017 Martin Sandve Alnæs
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Tools for precomputed tables of terminal values."""

import logging
import typing

import numpy

import ufl
import ufl.utils.derivativetuples
from ffcx.element_interface import basix_index, convert_element, QuadratureElement
from ffcx.ir.representationutils import (create_quadrature_points_and_weights,
                                         integral_type_to_entity_dim,
                                         map_integral_points)

logger = logging.getLogger("ffcx")

# Using same defaults as numpy.allclose
default_rtol = 1e-6
default_atol = 1e-9

piecewise_ttypes = ("piecewise", "fixed", "ones", "zeros")
uniform_ttypes = ("fixed", "ones", "zeros", "uniform")


class ModifiedTerminalElement(typing.NamedTuple):
    element: ufl.FiniteElementBase
    averaged: str
    local_derivatives: typing.Tuple[int]
    fc: int


class UniqueTableReferenceT(typing.NamedTuple):
    name: str
    values: numpy.typing.NDArray[numpy.float64]
    offset: int
    block_size: int
    ttype: str
    is_piecewise: bool
    is_uniform: bool
    is_permuted: bool


def equal_tables(a, b, rtol=default_rtol, atol=default_atol):
    a = numpy.asarray(a)
    b = numpy.asarray(b)
    if a.shape != b.shape:
        return False
    else:
        return numpy.allclose(a, b, rtol=rtol, atol=atol)


def clamp_table_small_numbers(table,
                              rtol=default_rtol,
                              atol=default_atol,
                              numbers=(-1.0, 0.0, 1.0)):
    """Clamp almost 0,1,-1 values to integers. Returns new table."""
    # Get shape of table and number of columns, defined as the last axis
    table = numpy.asarray(table)
    for n in numbers:
        table[numpy.where(numpy.isclose(table, n, rtol=rtol, atol=atol))] = n
    return table


def get_ffcx_table_values(points, cell, integral_type, element, avg, entitytype,
                          derivative_counts, flat_component):
    """Extract values from FFCx element table.

    Returns a 3D numpy array with axes
    (entity number, quadrature point number, dof number)
    """
    element = convert_element(element)
    deriv_order = sum(derivative_counts)

    if integral_type in ufl.custom_integral_types:
        # Use quadrature points on cell for analysis in custom integral types
        integral_type = "cell"
        assert not avg

    if integral_type == "expression":
        # FFCx tables for expression are generated as interior cell points
        integral_type = "cell"

    if avg in ("cell", "facet"):
        # Redefine points to compute average tables

        # Make sure this is not called with points, that doesn't make sense
        # assert points is None

        # Not expecting derivatives of averages
        assert not any(derivative_counts)
        assert deriv_order == 0

        # Doesn't matter if it's exterior or interior facet integral,
        # just need a valid integral type to create quadrature rule
        if avg == "cell":
            integral_type = "cell"
        elif avg == "facet":
            integral_type = "exterior_facet"

        if isinstance(element, QuadratureElement):
            points = element._points
            weights = element._weights
        else:
            # Make quadrature rule and get points and weights
            points, weights = create_quadrature_points_and_weights(
                integral_type, cell, element.highest_degree(), "default")

    # Tabulate table of basis functions and derivatives in points for each entity
    tdim = cell.topological_dimension()
    entity_dim = integral_type_to_entity_dim(integral_type, tdim)
    num_entities = ufl.cell.num_cell_entities[cell.cellname()][entity_dim]

    # Extract arrays for the right scalar component
    component_tables = []
    component_element, offset, stride = element.get_component_element(flat_component)

    for entity in range(num_entities):
        entity_points = map_integral_points(points, integral_type, cell, entity)
        tbl = component_element.tabulate(deriv_order, entity_points)
        tbl = tbl[basix_index(derivative_counts)]
        component_tables.append(tbl)

    if avg in ("cell", "facet"):
        # Compute numeric integral of the each component table
        wsum = sum(weights)
        for entity, tbl in enumerate(component_tables):
            num_dofs = tbl.shape[1]
            tbl = numpy.dot(tbl, weights) / wsum
            tbl = numpy.reshape(tbl, (1, num_dofs))
            component_tables[entity] = tbl

    # Loop over entities and fill table blockwise (each block = points x dofs)
    # Reorder axes as (points, dofs) instead of (dofs, points)
    assert len(component_tables) == num_entities
    num_points, num_dofs = component_tables[0].shape
    shape = (1, num_entities, num_points, num_dofs)
    res = numpy.zeros(shape)
    for entity in range(num_entities):
        res[:, entity, :, :] = component_tables[entity]

    return {'array': res, 'offset': offset, 'stride': stride}


def generate_psi_table_name(quadrature_rule, element_counter, averaged: str, entitytype, derivative_counts,
                            flat_component):
    """Generate a name for the psi table.

    Format:
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

    Q   - unique ID of quadrature rule, to distinguish between tables in a mixed quadrature rule setting

    """
    name = "FE%d" % element_counter
    if flat_component is not None:
        name += "_C%d" % flat_component
    if any(derivative_counts):
        name += "_D" + "".join(str(d) for d in derivative_counts)
    name += {None: "", "cell": "_AC", "facet": "_AF"}[averaged]
    name += {"cell": "", "facet": "_F", "vertex": "_V"}[entitytype]
    name += f"_Q{quadrature_rule.id()}"
    return name


def get_modified_terminal_element(mt) -> typing.Optional[ModifiedTerminalElement]:
    gd = mt.global_derivatives
    ld = mt.local_derivatives

    # Extract element from FormArguments and relevant GeometricQuantities
    if isinstance(mt.terminal, ufl.classes.FormArgument):
        if gd and mt.reference_value:
            raise RuntimeError(
                "Global derivatives of reference values not defined.")
        elif ld and not mt.reference_value:
            raise RuntimeError(
                "Local derivatives of global values not defined.")
        element = convert_element(mt.terminal.ufl_function_space().ufl_element())
        fc = mt.flat_component
    elif isinstance(mt.terminal, ufl.classes.SpatialCoordinate):
        if mt.reference_value:
            raise RuntimeError("Not expecting reference value of x.")
        if gd:
            raise RuntimeError("Not expecting global derivatives of x.")
        element = convert_element(mt.terminal.ufl_domain().ufl_coordinate_element())
        if not ld:
            fc = mt.flat_component
        else:
            # Actually the Jacobian expressed as reference_grad(x)
            fc = mt.flat_component  # x-component
            assert len(mt.component) == 1
            assert mt.component[0] == mt.flat_component
    elif isinstance(mt.terminal, ufl.classes.Jacobian):
        if mt.reference_value:
            raise RuntimeError("Not expecting reference value of J.")
        if gd:
            raise RuntimeError("Not expecting global derivatives of J.")
        element = convert_element(mt.terminal.ufl_domain().ufl_coordinate_element())
        assert len(mt.component) == 2
        # Translate component J[i,d] to x element context rgrad(x[i])[d]
        fc, d = mt.component  # x-component, derivative
        ld = tuple(sorted((d, ) + ld))
    else:
        return None

    assert (mt.averaged is None) or not (ld or gd)
    # Change derivatives format for table lookup
    gdim = mt.terminal.ufl_domain().geometric_dimension()
    local_derivatives = ufl.utils.derivativetuples.derivative_listing_to_counts(
        ld, gdim)

    return ModifiedTerminalElement(element, mt.averaged, local_derivatives, fc)


def permute_quadrature_interval(points, reflections=0):
    output = points.copy()
    for p in output:
        assert len(p) < 2 or numpy.isclose(p[1], 0)
        assert len(p) < 3 or numpy.isclose(p[2], 0)
    for i in range(reflections):
        for n, p in enumerate(output):
            output[n] = [1 - p[0]]
    return output


def permute_quadrature_triangle(points, reflections=0, rotations=0):
    output = points.copy()
    for p in output:
        assert len(p) < 3 or numpy.isclose(p[2], 0)
    for i in range(rotations):
        for n, p in enumerate(output):
            output[n] = [p[1], 1 - p[0] - p[1]]
    for i in range(reflections):
        for n, p in enumerate(output):
            output[n] = [p[1], p[0]]
    return output


def permute_quadrature_quadrilateral(points, reflections=0, rotations=0):
    output = points.copy()
    for p in output:
        assert len(p) < 3 or numpy.isclose(p[2], 0)
    for i in range(rotations):
        for n, p in enumerate(output):
            output[n] = [p[1], 1 - p[0]]
    for i in range(reflections):
        for n, p in enumerate(output):
            output[n] = [p[1], p[0]]
    return output


def build_optimized_tables(quadrature_rule, cell, integral_type, entitytype,
                           modified_terminals, existing_tables,
                           rtol=default_rtol, atol=default_atol):
    """Build the element tables needed for a list of modified terminals.

    Input:
      entitytype - str
      modified_terminals - ordered sequence of unique modified terminals
      FIXME: Document

    Output:
      mt_tables - dict(ModifiedTerminal: table data)
    """
    # Add to element tables
    analysis = {}
    for mt in modified_terminals:
        res = get_modified_terminal_element(mt)
        if res:
            analysis[mt] = res

    # Build element numbering using topological ordering so subelements
    # get priority
    all_elements = [res[0] for res in analysis.values()]
    unique_elements = ufl.algorithms.sort_elements(
        ufl.algorithms.analysis.extract_sub_elements(all_elements))
    element_numbers = {element: i for i, element in enumerate(unique_elements)}
    mt_tables = {}

    _existing_tables = existing_tables.copy()

    for mt in modified_terminals:
        res = analysis.get(mt)
        if not res:
            continue
        element, avg, local_derivatives, flat_component = res

        # Generate table and store table name with modified terminal

        # Build name for this particular table
        element_number = element_numbers[element]
        name = generate_psi_table_name(quadrature_rule, element_number, avg, entitytype,
                                       local_derivatives, flat_component)

        # FIXME - currently just recalculate the tables every time,
        # only reusing them if they match numerically.
        # It should be possible to reuse the cached tables by name, but
        # the dofmap offset may differ due to restriction.

        tdim = cell.topological_dimension()

        if integral_type == "interior_facet":
            if tdim == 1:
                t = get_ffcx_table_values(quadrature_rule.points, cell,
                                          integral_type, element, avg, entitytype,
                                          local_derivatives, flat_component)
            elif tdim == 2:
                new_table = []
                for ref in range(2):
                    new_table.append(get_ffcx_table_values(
                        permute_quadrature_interval(quadrature_rule.points, ref), cell,
                        integral_type, element, avg, entitytype, local_derivatives, flat_component))

                t = new_table[0]
                t['array'] = numpy.vstack([td['array'] for td in new_table])
            elif tdim == 3:
                cell_type = cell.cellname()
                if cell_type == "tetrahedron":
                    new_table = []
                    for rot in range(3):
                        for ref in range(2):
                            new_table.append(get_ffcx_table_values(
                                permute_quadrature_triangle(
                                    quadrature_rule.points, ref, rot),
                                cell, integral_type, element, avg, entitytype, local_derivatives,
                                flat_component))
                    t = new_table[0]
                    t['array'] = numpy.vstack([td['array'] for td in new_table])
                elif cell_type == "hexahedron":
                    new_table = []
                    for rot in range(4):
                        for ref in range(2):
                            new_table.append(get_ffcx_table_values(
                                permute_quadrature_quadrilateral(
                                    quadrature_rule.points, ref, rot),
                                cell, integral_type, element, avg, entitytype, local_derivatives, flat_component))
                    t = new_table[0]
                    t['array'] = numpy.vstack([td['array'] for td in new_table])
        else:
            t = get_ffcx_table_values(quadrature_rule.points, cell,
                                      integral_type, element, avg, entitytype,
                                      local_derivatives, flat_component)
        # Clean up table
        tbl = clamp_table_small_numbers(t['array'], rtol=rtol, atol=atol)
        tabletype = analyse_table_type(tbl)

        if tabletype in piecewise_ttypes:
            # Reduce table to dimension 1 along num_points axis in generated code
            tbl = tbl[:, :, :1, :]
        if tabletype in uniform_ttypes:
            # Reduce table to dimension 1 along num_entities axis in generated code
            tbl = tbl[:, :1, :, :]
        is_permuted = is_permuted_table(tbl)
        if not is_permuted:
            # Reduce table along num_perms axis
            tbl = tbl[:1, :, :, :]

        # Check for existing identical table
        new_table = True
        for table_name in _existing_tables:
            if equal_tables(tbl, _existing_tables[table_name]):
                name = table_name
                tbl = _existing_tables[name]
                new_table = False
                break

        if new_table:
            _existing_tables[name] = tbl

        cell_offset = 0
        element = convert_element(element)

        if mt.restriction == "-" and isinstance(mt.terminal, ufl.classes.FormArgument):
            # offset = 0 or number of element dofs, if restricted to "-"
            cell_offset = element.dim

        offset = cell_offset + t['offset']
        block_size = t['stride']

        # tables is just np.arrays, mt_tables hold metadata too
        mt_tables[mt] = UniqueTableReferenceT(
            name, tbl, offset, block_size, tabletype,
            tabletype in piecewise_ttypes, tabletype in uniform_ttypes, is_permuted)

    return mt_tables


def is_zeros_table(table, rtol=default_rtol, atol=default_atol):
    return (numpy.product(table.shape) == 0
            or numpy.allclose(table, numpy.zeros(table.shape), rtol=rtol, atol=atol))


def is_ones_table(table, rtol=default_rtol, atol=default_atol):
    return numpy.allclose(table, numpy.ones(table.shape), rtol=rtol, atol=atol)


def is_quadrature_table(table, rtol=default_rtol, atol=default_atol):
    _, num_entities, num_points, num_dofs = table.shape
    Id = numpy.eye(num_points)
    return (num_points == num_dofs and all(
        numpy.allclose(table[0, i, :, :], Id, rtol=rtol, atol=atol) for i in range(num_entities)))


def is_permuted_table(table, rtol=default_rtol, atol=default_atol):
    return not all(
        numpy.allclose(table[0, :, :, :],
                       table[i, :, :, :], rtol=rtol, atol=atol)
        for i in range(1, table.shape[0]))


def is_piecewise_table(table, rtol=default_rtol, atol=default_atol):
    return all(
        numpy.allclose(table[0, :, 0, :],
                       table[0, :, i, :], rtol=rtol, atol=atol)
        for i in range(1, table.shape[2]))


def is_uniform_table(table, rtol=default_rtol, atol=default_atol):
    return all(
        numpy.allclose(table[0, 0, :, :],
                       table[0, i, :, :], rtol=rtol, atol=atol)
        for i in range(1, table.shape[1]))


def analyse_table_type(table, rtol=default_rtol, atol=default_atol):
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
