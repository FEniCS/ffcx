# Copyright (C) 2013-2017 Martin Sandve Alnæs
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Tools for precomputed tables of terminal values."""

import collections
import logging

import numpy

import ufl
import ufl.utils.derivativetuples
from ffcx.libtab_interface import create_libtab_element, libtab_index
from ffcx.ir.representationutils import (create_quadrature_points_and_weights,
                                         integral_type_to_entity_dim,
                                         map_integral_points)

logger = logging.getLogger("ffcx")

# Using same defaults as numpy.allclose
default_rtol = 1e-5
default_atol = 1e-8

table_origin_t = collections.namedtuple(
    "table_origin", ["element", "avg", "derivatives", "flat_component", "dofrange", "dofmap"])

piecewise_ttypes = ("piecewise", "fixed", "ones", "zeros")
uniform_ttypes = ("fixed", "ones", "zeros", "uniform")

valid_ttypes = set(("quadrature", )) | set(
    piecewise_ttypes) | set(uniform_ttypes)

unique_table_reference_t = collections.namedtuple(
    "unique_table_reference",
    ["name", "values", "dofrange", "dofmap", "original_dim", "ttype", "is_piecewise", "is_uniform",
     "is_permuted", "dof_base_permutations", "needs_permutation_data"])


# TODO: Get restriction postfix from somewhere central
def ufc_restriction_offset(restriction, length):
    if restriction == "-":
        return length
    else:
        return 0


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
                              numbers=(-1.0, -0.5, 0.0, 0.5, 1.0)):
    """Clamp almost 0,1,-1 values to integers. Returns new table."""
    # Get shape of table and number of columns, defined as the last axis
    table = numpy.asarray(table)
    for n in numbers:
        table[numpy.where(numpy.isclose(table, n, rtol=rtol, atol=atol))] = n
    return table


def strip_table_zeros(table, block_size, rtol=default_rtol, atol=default_atol):
    """Strip zero columns from table. Returns column range (begin, end) and the new compact table."""
    # Get shape of table and number of columns, defined as the last axis
    table = numpy.asarray(table)
    sh = table.shape

    # Find nonzero columns
    z = numpy.zeros(sh[:-1])  # Correctly shaped zero table
    dofmap = tuple(
        i for i in range(sh[-1]) if not numpy.allclose(z, table[..., i], rtol=rtol, atol=atol))
    if dofmap:
        # Find first nonzero column
        begin = dofmap[0]
        # Find (one beyond) last nonzero column
        end = dofmap[-1] + 1
    else:
        begin = 0
        end = 0

    for i in dofmap:
        if i % block_size != dofmap[0] % block_size:
            # If dofs are not all in the same block component, don't remove intermediate zeros
            dofmap = tuple(range(begin, end))
            break
    else:
        # If dofs are all in the same block component, keep only that block component
        dofmap = tuple(range(begin, end, block_size))

    # Make subtable by dropping zero columns
    stripped_table = table[..., dofmap]
    dofrange = (begin, end)
    return dofrange, dofmap, stripped_table


def build_unique_tables(tables, rtol=default_rtol, atol=default_atol):
    """Return list of unique tables.

    Given a list or dict of tables, return a list of unique tables
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


def get_ffcx_table_values(points, cell, integral_type, ufl_element, avg, entitytype,
                          derivative_counts, flat_component):
    """Extract values from ffcx element table.

    Returns a 3D numpy array with axes
    (entity number, quadrature point number, dof number)
    """
    deriv_order = sum(derivative_counts)

    if integral_type in ufl.custom_integral_types:
        # Use quadrature points on cell for analysis in custom integral types
        integral_type = "cell"
        assert not avg

    if integral_type == "expression":
        # FFCX tables for expression are generated as interior cell points
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
    libtab_element = create_libtab_element(ufl_element)

    # Extract arrays for the right scalar component
    component_tables = []
    sh = ufl_element.value_shape()
    if sh == ():
        # Scalar valued element
        for entity in range(num_entities):
            entity_points = map_integral_points(
                points, integral_type, cell, entity)
            # libtab
            tbl = libtab_element.tabulate(deriv_order, entity_points)
            index = libtab_index(*derivative_counts)
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
            tbl = libtab_element.tabulate(deriv_order, entity_points)
            tbl = tbl[libtab_index(*derivative_counts)]
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
        sub_dims = [0] + [e.dim for e in libtab_element.sub_elements]
        sub_cmps = [0] + [e.value_size for e in libtab_element.sub_elements]

        irange = numpy.cumsum(sub_dims)
        crange = numpy.cumsum(sub_cmps)

        # Find index of sub element which corresponds to the current flat component
        component_element_index = numpy.where(
            crange <= flat_component)[0].shape[0] - 1

        ir = irange[component_element_index:component_element_index + 2]
        cr = crange[component_element_index:component_element_index + 2]

        component_element = libtab_element.sub_elements[component_element_index]

        # Get the block size to switch XXYYZZ ordering to XYZXYZ
        if isinstance(ufl_element, ufl.VectorElement) or isinstance(ufl_element, ufl.TensorElement):
            block_size = libtab_element.block_size
            ir = [ir[0] * block_size // irange[-1], irange[-1], block_size]

        def slice_size(r):
            if len(r) == 1:
                return r[0]
            if len(r) == 2:
                return r[1] - r[0]
            if len(r) == 3:
                return 1 + (r[1] - r[0] - 1) // r[2]

        # Follows from FIAT's MixedElement tabulation
        # Tabulating MixedElement in FIAT would result in tabulated subelements
        # padded with zeros
        for entity in range(num_entities):
            entity_points = map_integral_points(
                points, integral_type, cell, entity)

            # libtab
            tbl = component_element.tabulate(
                deriv_order, entity_points)
            index = libtab_index(*derivative_counts)
            tbl = tbl[index].transpose()

            # Prepare a padded table with zeros
            padded_shape = (libtab_element.dim,) + libtab_element.value_shape + (len(entity_points), )
            padded_tbl = numpy.zeros(padded_shape, dtype=tbl.dtype)

            tab = tbl.reshape(slice_size(ir), slice_size(cr), -1)
            padded_tbl[slice(*ir), slice(*cr)] = tab

            component_tables.append(padded_tbl[:, flat_component, :])

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
                         rtol=default_rtol,
                         atol=default_atol):
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

    # Build element numbering using topological ordering so subelements
    # get priority
    all_elements = [res[0] for res in analysis.values()]
    unique_elements = ufl.algorithms.sort_elements(
        ufl.algorithms.analysis.extract_sub_elements(all_elements))
    element_numbers = {element: i for i, element in enumerate(unique_elements)}

    def add_table(res):
        element, avg, local_derivatives, flat_component = res

        # Build name for this particular table
        element_number = element_numbers[element]
        name = generate_psi_table_name(quadrature_rule, element_number, avg, entitytype,
                                       local_derivatives, flat_component)

        if name not in tables:
            tdim = cell.topological_dimension()
            if entitytype == "facet":
                if tdim == 1:
                    tables[name] = numpy.array([
                        get_ffcx_table_values(quadrature_rule.points, cell,
                                              integral_type, element, avg, entitytype,
                                              local_derivatives, flat_component)])
                elif tdim == 2:
                    # Extract the values of the table from ffc table format
                    new_table = []
                    for ref in range(2):
                        new_table.append(get_ffcx_table_values(
                            permute_quadrature_interval(
                                quadrature_rule.points, ref),
                            cell, integral_type, element, avg, entitytype, local_derivatives, flat_component))

                    tables[name] = numpy.array(new_table)
                elif tdim == 3:
                    cell_type = cell.cellname()
                    if cell_type == "tetrahedron":
                        # Extract the values of the table from ffc table format
                        new_table = []
                        for rot in range(3):
                            for ref in range(2):
                                new_table.append(get_ffcx_table_values(
                                    permute_quadrature_triangle(
                                        quadrature_rule.points, ref, rot),
                                    cell, integral_type, element, avg, entitytype, local_derivatives, flat_component))

                        tables[name] = numpy.array(new_table)
                    elif cell_type == "hexahedron":
                        # Extract the values of the table from ffc table format
                        new_table = []
                        for rot in range(4):
                            for ref in range(2):
                                new_table.append(get_ffcx_table_values(
                                    permute_quadrature_quadrilateral(
                                        quadrature_rule.points, ref, rot),
                                    cell, integral_type, element, avg, entitytype, local_derivatives, flat_component))

                        tables[name] = numpy.array(new_table)
            else:
                # Extract the values of the table from ffc table format
                tables[name] = numpy.array([get_ffcx_table_values(quadrature_rule.points, cell,
                                                                  integral_type, element, avg, entitytype,
                                                                  local_derivatives, flat_component)])

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
        # for subelement in sort_elements(extract_sub_elements([element])):
        #     for fc in range(product(subelement.reference_value_shape())):
        #         subres = (subelement, avg, local_derivatives, fc)
        #         name_ignored = add_table(subres)

        # Generate table and store table name with modified terminal
        name = add_table(res)
        mt_table_names[mt] = name

    return tables, mt_table_names, table_origins


def optimize_element_tables(tables,
                            table_origins,
                            rtol=default_rtol,
                            atol=default_atol):
    """Optimize tables and make unique set.

    Steps taken:

      - clamp values that are very close to -1, 0, +1 to those values

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
    table_permuted = {}
    table_original_num_dofs = {}

    for name in used_names:
        tbl = tables[name]

        # Clamp to selected small numbers if close,
        # (-1.0, -0.5, 0.0, 0.5 and +1.0)
        # (i.e. 0.999999 -> 1.0 if within rtol/atol distance)
        tbl = clamp_table_small_numbers(tbl, rtol=rtol, atol=atol)

        # Store original dof dimension before compressing
        num_dofs = tbl.shape[3]
        ufl_element = table_origins[name][0]
        block_size = 1
        if isinstance(ufl_element, ufl.VectorElement) or isinstance(ufl_element, ufl.TensorElement):
            block_size = len(ufl_element.sub_elements())

        dofrange, dofmap, tbl = strip_table_zeros(
            tbl, block_size, rtol=rtol, atol=atol)

        compressed_tables[name] = tbl
        table_ranges[name] = dofrange
        table_dofmaps[name] = dofmap
        table_permuted[name] = is_permuted_table(tbl)
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
    table_unames = {
        name: unique_names[name_to_unique_index[name]] for name in name_to_unique_index}

    # Build mapping from unique table name to the table itself
    unique_tables = {}
    unique_table_origins = {}
    for ui, tbl in enumerate(unique_tables_list):
        uname = unique_names[ui]
        unique_tables[uname] = tbl
        unique_table_origins[uname] = table_origins[uname]

    return unique_tables, unique_table_origins, table_unames, table_ranges, table_dofmaps, table_permuted, \
        table_original_num_dofs


def is_zeros_table(table, rtol=default_rtol, atol=default_atol):
    return (ufl.utils.sequences.product(table.shape) == 0
            or numpy.allclose(table, numpy.zeros(table.shape), rtol=rtol, atol=atol))


def is_ones_table(table, rtol=default_rtol, atol=default_atol):
    return numpy.allclose(table, numpy.ones(table.shape), rtol=rtol, atol=atol)


def is_quadrature_table(table, rtol=default_rtol, atol=default_atol):
    num_perms, num_entities, num_points, num_dofs = table.shape
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


def analyse_table_type(table, rtol=default_rtol, atol=default_atol):
    num_perms, num_entities, num_points, num_dofs = table.shape
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


def is_uniform_table(table, rtol=default_rtol, atol=default_atol):
    return all(
        numpy.allclose(table[0, 0, :, :],
                       table[0, i, :, :], rtol=rtol, atol=atol)
        for i in range(1, table.shape[1]))


def analyse_table_types(unique_tables, rtol=default_rtol, atol=default_atol):
    return {
        uname: analyse_table_type(table, rtol=rtol, atol=atol)
        for uname, table in unique_tables.items()
    }


def build_optimized_tables(quadrature_rule,
                           cell,
                           integral_type,
                           entitytype,
                           modified_terminals,
                           existing_tables,
                           rtol=default_rtol,
                           atol=default_atol):

    # Build tables needed by all modified terminals
    tables, mt_table_names, table_origins = build_element_tables(
        quadrature_rule,
        cell,
        integral_type,
        entitytype,
        modified_terminals,
        rtol=rtol,
        atol=atol)

    # Optimize tables and get table name and dofrange for each modified terminal
    unique_tables, unique_table_origins, table_unames, table_ranges, table_dofmaps, table_permuted, \
        table_original_num_dofs = optimize_element_tables(
            tables, table_origins, rtol=rtol, atol=atol)

    # Get num_dofs for all tables before they can be deleted later
    unique_table_num_dofs = {uname: tbl.shape[-1]
                             for uname, tbl in unique_tables.items()}

    # Analyze tables for properties useful for optimization
    unique_table_ttypes = analyse_table_types(
        unique_tables, rtol=rtol, atol=atol)

    # Compress tables that are constant along num_entities or num_points
    for uname, tabletype in unique_table_ttypes.items():
        if tabletype in piecewise_ttypes:
            # Reduce table to dimension 1 along num_points axis in generated code
            unique_tables[uname] = unique_tables[uname][:, :, :1, :]
        if tabletype in uniform_ttypes:
            # Reduce table to dimension 1 along num_entities axis in generated code
            unique_tables[uname] = unique_tables[uname][:, :1, :, :]
        if not table_permuted[uname]:
            # Reduce table to dimenstion 2 along num_perms axis in generated code
            unique_tables[uname] = unique_tables[uname][:1, :, :, :]

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

    needs_permutation_data = False
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
        is_permuted = table_permuted[name]
        if is_permuted:
            needs_permutation_data = True

        # Map name -> uname
        uname = table_unames[name]

        # Map uname -> ename
        ename = name_map.get(uname, uname)

        # Some more metadata stored under the ename
        ttype = unique_table_ttypes[ename]

        offset = 0
        # Add offset to dofmap and dofrange for restricted terminals
        if mt.restriction and isinstance(mt.terminal, ufl.classes.FormArgument):
            # offset = 0 or number of dofs before table optimization
            offset = ufc_restriction_offset(mt.restriction, original_dim)
            (b, e) = dofrange
            dofrange = (b + offset, e + offset)
            dofmap = tuple(i + offset for i in dofmap)

        base_perms = [
            [[p[i - offset][j - offset] for j in dofmap] for i in dofmap]
            for p in create_libtab_element(table_origins[name][0]).base_permutations]

        needs_permutation_data = False
        for p in base_perms:
            if not numpy.allclose(p, numpy.identity(len(p))):
                needs_permutation_data = True

        # Store reference to unique table for this mt
        mt_unique_table_reference[mt] = unique_table_reference_t(
            ename, unique_tables[ename], dofrange, dofmap, original_dim, ttype,
            ttype in piecewise_ttypes, ttype in uniform_ttypes, is_permuted, base_perms,
            needs_permutation_data)

    return (unique_tables, unique_table_ttypes, unique_table_num_dofs,
            mt_unique_table_reference)
