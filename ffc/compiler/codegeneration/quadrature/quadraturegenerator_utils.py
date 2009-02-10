"Utility functions for quadrature representation"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2007-03-16 -- 2008-09-08"
__copyright__ = "Copyright (C) 2007-2008 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Python modules
import os
from sets import Set

# FFC language modules
from ffc.compiler.language.index import *
from ffc.compiler.language.restriction import *

# FFC tensor representation modules
from ffc.compiler.representation.tensor.multiindex import *

# FFC code generation modules
from ffc.compiler.codegeneration.common.codegenerator import *

from numpy import shape, array

def compute_macro_idims(monomial, idims, irank):
    "Compute macro dimensions in case of restricted basisfunctions"

    # copy dims
    macro_idims = [dim for dim in idims]
    for i in range(irank):
        for v in monomial.basisfunctions:
            if v.index.type == Index.PRIMARY and v.index.index == i:
                if v.restriction != None:
                    macro_idims[i] = idims[i] * 2
                    break
    return macro_idims


def generate_psi_name(tensor_number, psi_indices, vindex, aindices, bindices,\
                      format, qeindices = []):
    "Generate the name of a Psi table, but not the access to the matrix."

    # Prefetch formats to speed up code generation
    format_secondary_index  = format["secondary index"]
    index_names             = format["psi index names"]
    multi_index = []

    indices = ""
    for index in psi_indices:
        # We only need QE psi tables if they're primary indices
        if index in qeindices and not index.type == Index.PRIMARY:
            return None
        # If we have a vindex we move it in front of the name
        if index == vindex:
            indices = format_secondary_index(index_names[index.type](index.index)) + indices
        else:
            indices += format_secondary_index(index_names[index.type](index([], aindices, bindices, [])))

    name = format["psis"] + format_secondary_index("t%d" %tensor_number) + indices

    return name

def generate_psi_entry2(num_ips, tensor_number, aindices, bindices, psi_indices,\
                        vindices, name_map, optimise_level, format, qeindices = []):
    "Generate both name and matrix access for a given Psi."

    # Prefetch formats to speed up code generation
    primary_indices         = [format["first free index"], format["second free index"]]
    format_secondary_index  = format["secondary index"]
    format_ip               = format["integration points"]
    index_names             = format["psi index names"]

    # If we only have 1 IP, pick first
    # FIXME: this is language dependent, other languages might use '1' as first
    # component
    if num_ips == 1:
        format_ip = "0"

    entry = ""
    indices = ""
    dof_range = -1
    loop_dof = 0
    df_range = 0
    # Create indices (for name) and information to generate loop
    for index in psi_indices:
        if index in vindices:
            indices = format_secondary_index(index_names[index.type](index.index)) + indices
            loop_dof = index(primary_indices, aindices, bindices, [])
            df_range = len(index.range)
        else:
            indices += format_secondary_index(index_names[index.type](index([], aindices, bindices, [])))

    name = format["psis"] + format_secondary_index("t%d" %tensor_number) + indices

    # For QE, we only need a psi entry if we do have a primary index
    for index in psi_indices:
        if index in qeindices and not index.type == Index.PRIMARY:
            name = ""

    # If we have a name and a name map, m
    if name and name in name_map:
        dof_num = loop_dof
        if name_map[name][1]:
            # Get non-zero column info and new name
            i, cols = name_map[name][1]
            name = name_map[name][0]
            dof_range = len(cols)
            if dof_range == 1:
                # Direct lookup
                dof_num = "%d" % cols[0]
                entry = name + format["matrix access"](format_ip, dof_num)
            else:
                dof_num = format["nonzero columns"](i) + format["array access"](loop_dof)
                entry = name + format["matrix access"](format_ip, loop_dof)
            if not name:
                entry = name
        else:
            name = name_map[name][0]
            if name:
                entry = name + format["matrix access"](format_ip, loop_dof)
            else:
                entry = name
    else:
        dof_num = loop_dof
        entry = ""

    # If optimise level is 0, we loop the entire range of the dof
    if optimise_level == 0:
        dof_num = loop_dof
        dof_range = -1

    return ([loop_dof, dof_num], dof_range, entry)

def generate_loop(lines, loop_vars, Indent, format):
    "This function generates a loop over a vector or matrix."

    # Prefetch formats to speed up code generation
    format_loop     = format["loop"]
    format_begin    = format["block begin"]
    format_end      = format["block end"]
    format_comment  = format["comment"]

    if not loop_vars:
        return lines
    code = []
    for ls in loop_vars:
        # Get index and lower and upper bounds
        index, lower, upper = ls

        # Loop index
        code.append( Indent.indent(format_loop(index, lower, upper)) )
        code.append( Indent.indent(format_begin) )

        # Increase indentation
        Indent.increase()

        # If this is the last loop, write values
        if index == loop_vars[-1][0]:
            for l in lines:
                if (isinstance(l,tuple) or isinstance(l,list)) and len(l) == 2:
                    code.append((Indent.indent(l[0]), l[1]))
                elif isinstance(l,str):
                    code.append(Indent.indent(l))
                else:
                    print "line: ", l
                    raise RuntimeError, "Line must be a string or a list or tuple of length 2"

    # Decrease indentation and write end blocks
    vrs = [lv[0] for lv in loop_vars]
    vrs.reverse()
    for lv in vrs:
        Indent.decrease()
        code.append( Indent.indent(format_end + format_comment("end loop over '%s'" %lv) ) )

    return code

def equal_loops(tensors):
    "Group tensor with an equal number of quadrature points and primary indices"

    # Loop all tensors and group number of quadrature points
    group_tensors = {}
    for i in range(len(tensors)):
        tens = tensors[i]
        num_points = len(tens.quadrature.weights)
        if num_points in group_tensors:
            group_tensors[num_points].append(i)
        else:
            group_tensors[num_points] = [i]

    return group_tensors

def get_names_tables(tensor, tensor_number, format):
    "Tabulate values of basis functions and their derivatives at quadrature points"

    tables = {}

    # Loop psis
    for psi in tensor.Psis:

        # Get values of psi, list of indices and index of basisfunction for psi
        values, indices, vindex = psi[0], psi[1], psi[2]

        # Create list of secondary indices
        sec_indices = [index for index in indices if not index == vindex]

        # Get the number of dofs and quadrature points
        num_dofs = len(vindex.range)
        num_quadrature_points = len(tensor.quadrature.weights)

        # Get lists of secondary indices and auxiliary_0 indices (global)
        aindices = tensor.a.indices
        b0indices = tensor.b0.indices

        # Reduce multi indices a and b0 (number of pertubations)
        # according to current psi
        a_reduced = [[0],]*len(aindices[0])
        for index in indices:
            if index.type == Index.SECONDARY and index != vindex:
                a_reduced[index.index] = range(len(index.range))
        a_reduced = MultiIndex(a_reduced).indices

        b0_reduced = [[0],]*len(b0indices[0])
        for index in indices:
            if index.type == Index.AUXILIARY_0 and index != vindex:
                b0_reduced[index.index] = range(len(index.range))
        b0_reduced = MultiIndex(b0_reduced).indices

        # Loop secondary and auxiliary indices to generate names and entries
        # in the psi table
        for a in a_reduced:
            for b in b0_reduced:
                # Generate name
                name = generate_psi_name(tensor_number, indices, vindex,\
                                      a, b, format, tensor.qei)
                if name:
                    multi_index = []
                    for index in sec_indices:
                        multi_index.append(index([], a, b, []))

                    # Get values from psi tensor, should have format values[dofs][quad_points]
                    vals = values[tuple(multi_index)]

                    # Check if the values have the correct dimensions, otherwise quadrature is not correctly
                    # implemented for the given form!!
                    if numpy.shape(vals) != (num_dofs, num_quadrature_points):
                        raise RuntimeError, "Quadrature is not correctly implemented for the given form!"

                    # Generate values (FIAT returns [dof, quad_points] transpose to [quad_points, dof])
                    value = numpy.transpose(vals)
                    tables[name] = value

    return tables

def unique_tables(tables, format_epsilon):
    """Removes tables with redundant values and returns a name_map and a
    inverse_name_map. E.g.,

    tables = {a:[0,1,2], b:[0,2,3], c:[0,1,2], d:[0,1,2]}
    results in:
    tables = {a:[0,1,2], b:[0,2,3]}
    name_map = {a:[c,d]}
    inverse_name_map = {a:a, b:b, c:a, d:a}"""

    name_map = {}
    inverse_name_map = {}
    names = [name for name in tables]
    mapped = []

    # Loop all tables to see if some are redundant
    for i in range(len(names)):
        name0 = names[i]
        if name0 in mapped:
            continue
        val0 = array(tables[name0])

        for j in range(i+1, len(names)):
            name1 = names[j]
            if name1 in mapped:
                continue
            val1 = array(tables[name1])

            # Check if dimensions match
            if shape(val0) == shape(val1):
                # Check if values are the same
                if sqrt(((val0 - val1)*(val0 - val1)).sum()) < format_epsilon:
                    mapped.append(name1)
                    del tables[name1]
                    if name0 in name_map:
                        name_map[name0].append(name1)
                    else:
                        name_map[name0] = [name1]
                    # Create inverse name map
                    inverse_name_map[name1] = name0

    # Add self
    for name in tables:
        if not name in inverse_name_map:
            inverse_name_map[name] = name

    return (name_map, inverse_name_map)

def get_ones(tables, format_epsilon):
    "Return names of tables for which all values are 1.0"
    names = []
    for name in tables:
        vals = tables[name]
        one = True
        for r in range(numpy.shape(vals)[0]):
            for c in range(numpy.shape(vals)[1]):
                if abs(vals[r][c] - 1.0) > format_epsilon:
                    one = False
        if one:
            names.append(name)
    return names

def contains_zeros(tables, format_epsilon):
    "Checks if any tables contains all zeros"

    for name in tables:
        vals = tables[name]
        zero = True
        for r in range(numpy.shape(vals)[0]):
            if not zero:
                break
            for c in range(numpy.shape(vals)[1]):
                # If just one value is different from zero, break loops
                if abs(vals[r][c]) > format_epsilon:
                    zero = False
                    break

        if zero:
            debug("\n*** Warning: this table only contains zeros. This is not critical,")
            debug("but it might slow down the runtime performance of your code!")
            debug("Do you take derivatives of a constant?\n")

def unique_psi_tables(tables, optimisation_level, format):
    "Determine if some tensors have the same tables (and same names)"

    # Get formats
    format_epsilon = format["epsilon"]

    # Get unique tables
    name_map, inverse_name_map = unique_tables(tables, format_epsilon)

    # Exclude tables with all ones
    names = get_ones(tables, format_epsilon)
    for name in names:
        del tables[name]
        if name in name_map:
            maps = name_map[name]
            for m in maps:
                del inverse_name_map[m]
        del inverse_name_map[name]

    # Set values to zero if they are lower than threshold
    for name in tables:
        # Get values
        vals = tables[name]
        for r in range(numpy.shape(vals)[0]):
            for c in range(numpy.shape(vals)[1]):
                if abs(vals[r][c]) < format_epsilon:
                    vals[r][c] = 0
        tables[name] = vals

    # Extract the column numbers that are non-zero
    # (only for optimisations higher than 0)
    i = 0
    non_zero_columns = {}
    if optimisation_level > 0:
        for name in tables:
            # Get values
            vals = tables[name]

            # Use the first row as reference
            non_zeros = list(vals[0].nonzero()[0])

            # If all columns in the first row are non zero, there's no point
            # in continuing
            if len(non_zeros) == numpy.shape(vals)[1]:
                continue

            # If we only have one row (IP) we just need the nonzero columns
            if numpy.shape(vals)[0] == 1:
                if list(non_zeros):
                    non_zeros.sort()
                    non_zero_columns[name] = (i, non_zeros)

                    # Possibly compress values
                    if optimisation_level == 0:
                        tables[name] = vals
                    else:
                        tables[name] = vals[:, non_zeros]
                    i += 1

            # Check if the remaining rows are nonzero in the same positions, else expand
            else:
                for j in range(numpy.shape(vals)[0] - 1):
                    # All rows must have the same non-zero columns
                    # for the optimization to work (at this stage)
                    new_non_zeros = list(vals[j+1].nonzero()[0])
                    if non_zeros != new_non_zeros:
                        non_zeros = non_zeros + [c for c in new_non_zeros if not c in non_zeros]
                        # If this results in all columns being non-zero, continue.
                        if len(non_zeros) == numpy.shape(vals)[1]:
                            continue

                # Only add nonzeros if it implies a reduction of columns
                if len(non_zeros) != numpy.shape(vals)[1]:
                    non_zeros.sort()
                    non_zero_columns[name] = (i, non_zeros)

                    # Compress values
                    tables[name] = vals[:, non_zeros]
                    i += 1

    # Add non-zero column info to inverse_name_map
    # (so we only need to pass around one name_map to code generating functions)
    for name in inverse_name_map:
        if inverse_name_map[name] in non_zero_columns:
            nzc = non_zero_columns[inverse_name_map[name]]
            inverse_name_map[name] = [inverse_name_map[name], nzc]
        else:
            inverse_name_map[name] = [inverse_name_map[name], ()]

    # Check if we have some zeros in the tables
    # FIXME: this shouldn't be necessary once UFL has removed all zero terms
    contains_zeros(tables, format_epsilon)

    # If we found non zero columns we might be able to reduce number of tables
    # further if optimise level is higher than 0
    if non_zero_columns:

        # Try reducing the tables. This is possible if some tables have become
        # identical as a consequence of compressing the tables
        nm, inv_nm = unique_tables(tables, format_epsilon)

        # Update name maps
        for name in inverse_name_map:
            if inverse_name_map[name][0] in inv_nm:
                inverse_name_map[name][0] = inv_nm[inverse_name_map[name][0]]
        for name in nm:
            maps = nm[name]
            for m in maps:
                if not name in name_map:
                    name_map[name] = []
                if m in name_map:
                    name_map[name] += name_map[m] + [m]
                    del name_map[m]
                else:
                    name_map[name].append(m)

        # Exclude tables with all ones
        names = get_ones(tables, format_epsilon)
        # Because these tables now contain ones as a consequence of compression
        # we still need to consider the non-zero columns when looking up values
        # in coefficient arrays. The psi entries can however we neglected and we
        # don't need to tabulate the values
        for name in names:
            if name in name_map:
                maps = name_map[name]
                for m in maps:
                    inverse_name_map[m][0] = ""
            if name in inverse_name_map:
                inverse_name_map[name][0] = ""
            tables[name] = None

    # Write protect info
    for name in inverse_name_map:
        inverse_name_map[name] = tuple(inverse_name_map[name])

    return (inverse_name_map, tables)

def unique_weight_tables(tensors, format):
    "Determine if some tensors have the same tables (and same names)"

    tables = {}
    # Loop tensors and get all tables
    for tensor_number in range(len(tensors)):
        weights = tensors[tensor_number].quadrature.weights
        tables[tensor_number] = weights

    name_map, inverse_name_map = unique_tables(tables, format["epsilon"])

    return (name_map, tables)

