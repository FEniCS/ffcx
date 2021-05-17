# Copyright (C) 2013-2017 Martin Sandve Alnæs
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Tools for precomputed tables of terminal values."""

import collections
import logging

import numpy

import ufl
import ufl.utils.derivativetuples
from ffcx.element_interface import create_element, basix_index
from ffcx.ir.representationutils import (create_quadrature_points_and_weights,
                                         integral_type_to_entity_dim,
                                         map_integral_points)

logger = logging.getLogger("ffcx")

# Using same defaults as numpy.allclose
default_rtol = 1e-6
default_atol = 1e-9

piecewise_ttypes = ("piecewise", "fixed", "ones", "zeros")
uniform_ttypes = ("fixed", "ones", "zeros", "uniform")

unique_table_reference_t = collections.namedtuple(
    "unique_table_reference",
    ["name", "values", "dofrange", "dofmap", "ttype", "is_piecewise", "is_uniform",
     "is_permuted", "dof_base_transformations", "needs_transformation_data"])


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


def get_ffcx_table_values(points, cell, integral_type, ufl_element, avg, entitytype,
                          derivative_counts, flat_component):
    """Extract values from FFCx element table.

    Returns a 3D numpy array with axes
    (entity number, quadrature point number, dof number)
    """
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

        # Make quadrature rule and get points and weights
        points, weights = create_quadrature_points_and_weights(integral_type, cell,
                                                               ufl_element.degree(), "default")

    # Tabulate table of basis functions and derivatives in points for each entity
    tdim = cell.topological_dimension()
    entity_dim = integral_type_to_entity_dim(integral_type, tdim)
    num_entities = ufl.cell.num_cell_entities[cell.cellname()][entity_dim]

    numpy.set_printoptions(suppress=True, precision=2)
    basix_element = create_element(ufl_element)

    offset = 0
    stride = 1

    # Extract arrays for the right scalar component
    component_tables = []
    sh = ufl_element.value_shape()
    if sh == ():
        # Scalar valued element
        for entity in range(num_entities):
            entity_points = map_integral_points(
                points, integral_type, cell, entity)
            # basix
            tbl = basix_element.tabulate(deriv_order, entity_points)
            index = basix_index(*derivative_counts)
            tbl = tbl[index].transpose()

            component_tables.append(tbl)
    elif len(sh) > 0 and ufl_element.num_sub_elements() == 0:
        # 2-tensor-valued elements, not a tensor product
        # mapping flat_component back to tensor component
        (_, f2t) = ufl.permutation.build_component_numbering(
            sh, ufl_element.symmetry())
        t_comp = f2t[flat_component]

        for entity in range(num_entities):
            entity_points = map_integral_points(
                points, integral_type, cell, entity)
            tbl = basix_element.tabulate(deriv_order, entity_points)
            tbl = tbl[basix_index(*derivative_counts)]
            sum_sh = sum(sh)
            bshape = (tbl.shape[0],) + sh + (tbl.shape[1] // sum_sh,)
            tbl = tbl.reshape(bshape).transpose()

            if len(sh) == 1:
                component_tables.append(tbl[:, t_comp[0], :])
            elif len(sh) == 2:
                component_tables.append(tbl[:, t_comp[0], t_comp[1], :])
            else:
                raise RuntimeError(
                    "Cannot tabulate tensor valued element with rank > 2")
    else:
        # Vector-valued or mixed element
        sub_dims = [0] + [e.dim for e in basix_element.sub_elements]
        sub_cmps = [0] + [e.value_size for e in basix_element.sub_elements]

        irange = numpy.cumsum(sub_dims)
        crange = numpy.cumsum(sub_cmps)

        # Find index of sub element which corresponds to the current flat component
        component_element_index = numpy.where(
            crange <= flat_component)[0].shape[0] - 1

        ir = irange[component_element_index:component_element_index + 2]
        cr = crange[component_element_index:component_element_index + 2]

        component_element = basix_element.sub_elements[component_element_index]

        # Get the block size to switch XXYYZZ ordering to XYZXYZ
        if isinstance(ufl_element, ufl.VectorElement) or isinstance(ufl_element, ufl.TensorElement):
            block_size = basix_element.block_size
            ir = [ir[0] * block_size // irange[-1], irange[-1], block_size]
            stride = block_size

        offset = ir[0]

        def slice_size(r):
            if len(r) == 1:
                return r[0]
            if len(r) == 2:
                return r[1] - r[0]
            if len(r) == 3:
                return 1 + (r[1] - r[0] - 1) // r[2]

        for entity in range(num_entities):
            entity_points = map_integral_points(
                points, integral_type, cell, entity)

            # basix
            tbl = component_element.tabulate(
                deriv_order, entity_points)
            index = basix_index(*derivative_counts)
            tbl = tbl[index].transpose()

            # Prepare a padded table with zeros
            # padded_shape = (slice_size(ir),) + basix_element.value_shape + (len(entity_points), )
            # padded_tbl = numpy.zeros(padded_shape, dtype=tbl.dtype)

            if basix_element.family_name == "mixed element" and not component_element.is_blocked:
                tab = numpy.zeros_like(tbl.reshape(slice_size(ir), slice_size(cr), -1))
                for basis_i in range(tab.shape[0]):
                    for value_i in range(tab.shape[1]):
                        tab[basis_i, value_i] = tbl[value_i * tab.shape[0] + basis_i]
            else:
                tab = tbl.reshape(slice_size(ir), slice_size(cr), -1)

            #  print(padded_shape, tab.shape, slice(*cr), flat_component)
            #   padded_tbl[:, slice(*cr), :] = tab
            c = flat_component - cr[0]
            # print(flat_component, 'offset=', offset, 'stride=', stride)
            component_tables.append(tab[:, c, :])

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
    shape = (1, num_entities, num_points, num_dofs)
    res = numpy.zeros(shape)
    for entity in range(num_entities):
        res[:, entity, :, :] = numpy.transpose(component_tables[entity])

    return {'array': res, 'offset': offset, 'stride': stride}


def generate_psi_table_name(quadrature_rule, element_counter, averaged, entitytype, derivative_counts,
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


def get_modified_terminal_element(mt):
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
        element = mt.terminal.ufl_element()
        fc = mt.flat_component
    elif isinstance(mt.terminal, ufl.classes.SpatialCoordinate):
        if mt.reference_value:
            raise RuntimeError("Not expecting reference value of x.")
        if gd:
            raise RuntimeError("Not expecting global derivatives of x.")
        element = mt.terminal.ufl_domain().ufl_coordinate_element()
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

        element = mt.terminal.ufl_domain().ufl_coordinate_element()
        assert len(mt.component) == 2
        # Translate component J[i,d] to x element context rgrad(x[i])[d]
        fc, d = mt.component  # x-component, derivative

        # Grad(Jacobian(...)) should be a local derivative
        ld = tuple(sorted((d, ) + gd + ld))
    else:
        return None

    assert not (mt.averaged and (ld or gd))

    # Change derivatives format for table lookup
    # gdim = mt.terminal.ufl_domain().geometric_dimension()
    # global_derivatives = derivative_listing_to_counts(gd, gdim)

    # Change derivatives format for table lookup
    tdim = mt.terminal.ufl_domain().topological_dimension()
    local_derivatives = ufl.utils.derivativetuples.derivative_listing_to_counts(
        ld, tdim)

    return element, mt.averaged, local_derivatives, fc


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


def build_element_tables(quadrature_rule,
                         cell,
                         integral_type,
                         entitytype,
                         modified_terminals,
                         existing_tables,
                         rtol=default_rtol,
                         atol=default_atol):
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
        # FIXME: Use a namedtuple for res
        res = get_modified_terminal_element(mt)
        if res:
            analysis[mt] = res

    # Build element numbering using topological ordering so subelements
    # get priority
    all_elements = [res[0] for res in analysis.values()]
    unique_elements = ufl.algorithms.sort_elements(
        ufl.algorithms.analysis.extract_sub_elements(all_elements))
    element_numbers = {element: i for i, element in enumerate(unique_elements)}

    tables = existing_tables
    mt_tables = {}

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

        # FIXME - currently skip this, and just recalculate the tables every time,
        # only reusing them if they match numerically. mt_tables can only be reused if the
        # metadata (e.g. dofmap etc.) is also known to be the same.
        if False:
            # Find existing entry in mt_tables
            for k in mt_tables:
                if mt_tables[k]['name'] == name:
                    mt_tables[mt] = mt_tables[k]
                    break
        else:
            tdim = cell.topological_dimension()
            if entitytype == "facet":
                if tdim == 1:
                    t = get_ffcx_table_values(quadrature_rule.points, cell,
                                              integral_type, element, avg, entitytype,
                                              local_derivatives, flat_component)
                elif tdim == 2:
                    new_table = []
                    for ref in range(2):
                        new_table.append(get_ffcx_table_values(
                            permute_quadrature_interval(
                                quadrature_rule.points, ref),
                            cell, integral_type, element, avg, entitytype, local_derivatives, flat_component))
                    table_values = numpy.array([td['array'][0] for td in new_table])
                    t = new_table[0]
                    t['array'] = table_values
                elif tdim == 3:
                    cell_type = cell.cellname()
                    if cell_type == "tetrahedron":
                        new_table = []
                        for rot in range(3):
                            for ref in range(2):
                                new_table.append(get_ffcx_table_values(
                                    permute_quadrature_triangle(
                                        quadrature_rule.points, ref, rot),
                                    cell, integral_type, element, avg, entitytype, local_derivatives, flat_component))

                        table_values = numpy.array([td['array'][0] for td in new_table])
                        t = new_table[0]
                        t['array'] = table_values
                    elif cell_type == "hexahedron":
                        new_table = []
                        for rot in range(4):
                            for ref in range(2):
                                new_table.append(get_ffcx_table_values(
                                    permute_quadrature_quadrilateral(
                                        quadrature_rule.points, ref, rot),
                                    cell, integral_type, element, avg, entitytype, local_derivatives, flat_component))

                        table_values = numpy.array([td['array'][0] for td in new_table])
                        t = new_table[0]
                        t['array'] = table_values
            else:
                t = get_ffcx_table_values(quadrature_rule.points, cell,
                                          integral_type, element, avg, entitytype,
                                          local_derivatives, flat_component)

            # Clean up table
            tbl = clamp_table_small_numbers(t['array'], rtol=rtol, atol=atol)
            tabletype = analyse_table_type(tbl)
            t['ttype'] = tabletype
            if tabletype in piecewise_ttypes:
                # Reduce table to dimension 1 along num_points axis in generated code
                tbl = tbl[:, :, :1, :]
            if tabletype in uniform_ttypes:
                # Reduce table to dimension 1 along num_entities axis in generated code
                tbl = tbl[:, :1, :, :]

            t['is_permuted'] = is_permuted_table(tbl)
            if not t['is_permuted']:
                # Reduce table along num_perms axis
                tbl = tbl[:1, :, :, :]

            # Check for existing identical table
            xname_found = False
            for xname in tables:
                if equal_tables(tbl, tables[xname]):
                    xname_found = True
                    break

            if xname_found:
                # print('found existing table equiv to ', name, ' at ', xname)
                name = xname
                # Retrieve existing table
                t['array'] = tables[name]
            else:
                # Store new table
                t['array'] = tbl
                tables[name] = t['array']

            cell_offset = 0
            basix_element = create_element(element)

            if mt.restriction == "-" and isinstance(mt.terminal, ufl.classes.FormArgument):
                # offset = 0 or number of element dofs, if restricted to "-"
                cell_offset = basix_element.dim

            t['name'] = name
            t['is_piecewise'] = t['ttype'] in piecewise_ttypes
            t['is_uniform'] = t['ttype'] in uniform_ttypes
            num_dofs = t['array'].shape[3]
            dofmap = tuple(cell_offset + t['offset'] + i * t['stride'] for i in range(num_dofs))
            t['dofmap'] = dofmap

            t['base_transformations'] = [[[p[i - cell_offset][j - cell_offset] for j in dofmap]
                                         for i in dofmap]
                                         for p in basix_element.base_transformations]

            # tables is just np.arrays, mt_tables hold metadata too
            mt_tables[mt] = t

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


def build_optimized_tables(quadrature_rule,
                           cell,
                           integral_type,
                           entitytype,
                           modified_terminals,
                           existing_tables,
                           rtol=default_rtol,
                           atol=default_atol):

    # Build tables needed by all modified terminals
    mt_tables = build_element_tables(
        quadrature_rule,
        cell,
        integral_type,
        entitytype,
        modified_terminals,
        existing_tables,
        rtol=rtol,
        atol=atol)

    mt_unique_table_reference = {}
    for mt, table_data in mt_tables.items():

        # FIXME: Use offset and stride instead of dofmap and dofrange
        dofmap = table_data['dofmap']
        dofrange = (dofmap[0], dofmap[-1] + 1)
        is_permuted = table_data['is_permuted']

        needs_transformation_data = False
        # if is_permuted:
        #    needs_transformation_data = True

        base_transformations = table_data['base_transformations']
        for p in base_transformations:
            if not numpy.allclose(p, numpy.identity(len(p))):
                needs_transformation_data = True

        ttype = table_data['ttype']
        # Store reference to unique table for this mt
        mt_unique_table_reference[mt] = unique_table_reference_t(
            table_data['name'], table_data['array'], dofrange, dofmap, ttype,
            ttype in piecewise_ttypes, ttype in uniform_ttypes, is_permuted, base_transformations,
            needs_transformation_data)

    return mt_unique_table_reference
