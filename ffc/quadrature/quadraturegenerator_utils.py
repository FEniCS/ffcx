"Utility functions for quadrature representation."

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@gmail.com)"
__date__ = "2007-03-16"
__copyright__ = "Copyright (C) 2007-2010 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-01-21

# Python modules.
from numpy import transpose
from numpy import sqrt
from numpy import shape
from numpy import array

# FFC modules.
from ffc.log import debug
from ffc.log import error
from ffc.log import ffc_assert

def generate_loop(lines, loop_vars, Indent, format):
    "This function generates a loop over a vector or matrix."

    # Prefetch formats to speed up code generation.
    format_loop     = format["loop"]
    format_begin    = format["block begin"]
    format_end      = format["block end"]
    format_comment  = format["comment"]

    if not loop_vars:
        return lines
    code = []
    for ls in loop_vars:
        # Get index and lower and upper bounds.
        index, lower, upper = ls

        # Loop index.
        code.append( Indent.indent(format_loop(index, lower, upper)) )
        code.append( Indent.indent(format_begin) )

        # Increase indentation.
        Indent.increase()

        # If this is the last loop, write values.
        if index == loop_vars[-1][0]:
            for l in lines:
                if (isinstance(l,tuple) or isinstance(l,list)) and len(l) == 2:
                    code.append((Indent.indent(l[0]), l[1]))
                elif isinstance(l,str):
                    code.append(Indent.indent(l))
                else:
                    print "line: ", l
                    error("Line must be a string or a list or tuple of length 2.")

    # Decrease indentation and write end blocks.
    vrs = [lv[0] for lv in loop_vars]
    vrs.reverse()
    for lv in vrs:
        Indent.decrease()
        code.append( Indent.indent(format_end + format_comment("end loop over '%s'" %lv) ) )

    return code

def generate_psi_name(counter, facet, component, derivatives):
    """Generate a name for the psi table of the form:
    FE#_f#_C#_D###, where '#' will be an integer value.

    FE  - is a simple counter to distinguish the various bases, it will be
          assigned in an arbitrary fashion.

    f   - denotes facets if applicable, range(element.num_facets()).

    C   - is the component number if any (this does not yet take into account
          tensor valued functions)

    D   - is the number of derivatives in each spatial direction if any. If the
          element is defined in 3D, then D012 means d^3(*)/dydz^2."""

    name = "FE%d" % counter
    if not facet is None:
        name += "_f%d" % facet
    if component != () and component != []:
        name += "_C%d" % component
    if any(derivatives):
        name += "_D" + "".join([str(d) for d in derivatives])

    return name

def create_psi_tables(tables, format_epsilon, options):
    "Create names and maps for tables and non-zero entries if appropriate."

    debug("\nQG-utils, psi_tables:\n" + str(tables))
    # Create element map {points:{element:number,},}
    # and a plain dictionary {name:values,}.
    element_map, flat_tables = flatten_psi_tables(tables)
    debug("\nQG-utils, psi_tables, flat_tables:\n" + str(flat_tables))

    # Reduce tables such that we only have those tables left with unique values
    # Create a name map for those tables that are redundant.
    name_map, unique_tables = unique_psi_tables(flat_tables, format_epsilon, options)

    debug("\nQG-utils, psi_tables, unique_tables:\n" + str(unique_tables))
    debug("\nQG-utils, psi_tables, name_map:\n" + str(name_map))

    return (element_map, name_map, unique_tables)

def flatten_psi_tables(tables):
    """Create a 'flat' dictionary of tables with unique names and a name
    map that maps number of quadrature points and element name to a unique
    element number. returns:
    name_map    - {num_quad_points:{ufl_element:element_number,},}.
    flat_tables - {unique_table_name:values (ip,basis),}."""

    # Initialise return values and element counter.
    flat_tables = {}
    element_map = {}
    counter = 0
    # Loop quadrature points and get element dictionary {elem: {tables}}.
    for point in sorted(tables.keys()):
        elem_dict = tables[point]
        element_map[point] = {}
        debug("\nQG-utils, flatten_tables, points:\n" + str(point))
        debug("\nQG-utils, flatten_tables, elem_dict:\n" + str(elem_dict))
        # Loop all elements and get all their tables.
        for elem in sorted(elem_dict.keys(), lambda x, y: cmp(str(x), str(y))):
            facet_tables = elem_dict[elem]
            debug("\nQG-utils, flatten_tables, elem:\n" + str(elem))
            debug("\nQG-utils, flatten_tables, facet_tables:\n" + str(facet_tables))
            element_map[point][elem] = counter
            for facet in sorted(facet_tables.keys()):
                elem_table = facet_tables[facet]
                # If the element value rank != 0, we must loop the components.
                # before the derivatives (that's the way the values are tabulated).
                if len(elem.value_shape()) != 0:
                    for num_comp, comp_table in enumerate(elem_table):
                        for derivs in sorted(comp_table.keys()):
                            psi_table = num_deriv[derivs]
                            debug("\nQG-utils, flatten_tables, derivs:\n" + str(derivs))
                            debug("\nQG-utils, flatten_tables, psi_table:\n" + str(psi_table))
                            # Verify shape of basis (can be omitted for speed
                            # if needed I think).
                            ffc_assert(len(shape(psi_table)) == 2 and shape(psi_table)[1] == point, \
                                        "Something is wrong with this table: " + str(psi_table))
                            # Generate the table name.
                            name = generate_psi_name(counter, facet, num_comp, derivs)
                            debug("Table name: " + name)
                            ffc_assert(name not in flat_tables, \
                                        "Table name is not unique, something is wrong: " + name + str(flat_tables))
                            # Take transpose such that we get (ip_number, basis_number)
                            # instead of (basis_number, ip_number).
                            flat_tables[name] = transpose(psi_table)
                # If we don't have any components.
                else:
                    for derivs in sorted(elem_table.keys()):
                        psi_table = elem_table[derivs]
                        debug("\nQG-utils, flatten_tables, derivs:\n" + str(derivs))
                        debug("\nQG-utils, flatten_tables, psi_table:\n" + str(psi_table))
                        # Verify shape of basis (can be omitted for speed
                        # if needed I think).
                        ffc_assert(len(shape(psi_table)) == 2 and shape(psi_table)[1] == point, \
                                    "Something is wrong with this table: " + str(psi_table))
                        # Generate the table name.
                        name = generate_psi_name(counter, facet, (), derivs)
                        debug("Table name: " + name)
                        ffc_assert(name not in flat_tables, \
                                    "Table name is not unique, something is wrong: " + name + str(flat_tables))
                        flat_tables[name] = transpose(psi_table)
            # Increase unique element counter.
            counter += 1

    return (element_map, flat_tables)

def unique_psi_tables(tables, format_epsilon, options):
    """Returns a name map and a dictionary of unique tables. The function checks
    if values in the tables are equal, if this is the case it creates a name
    mapping. It also create additional information (depending on which options
    are set) such as if the table contains all ones, or only zeros, and a list
    on non-zero columns.
    unique_tables - {name:values,}.
    name_map      - {original_name:[new_name, non-zero-columns (list), is zero (bool), is ones (bool)],}."""

    # Get unique tables (from old table utility).
    name_map, inverse_name_map = unique_tables(tables, format_epsilon)

    debug("\ntables: " + str(tables))
    debug("\nname_map: " + str(name_map))
    debug("\ninv_name_map: " + str(inverse_name_map))

    # Get names of tables with all ones (from old table utility).
    names_ones = get_ones(tables, format_epsilon)

    # Set values to zero if they are lower than threshold.
    for name in tables:
        # Get values.
        vals = tables[name]
        for r in range(shape(vals)[0]):
            for c in range(shape(vals)[1]):
                if abs(vals[r][c]) < format_epsilon:
                    vals[r][c] = 0
        tables[name] = vals

    # Extract the column numbers that are non-zero.
    # If optimisation option is set
    # counter for non-zero column arrays.
    i = 0
    non_zero_columns = {}
    if options["non zero columns"]:
        for name in sorted(tables.keys()):
            # Get values.
            vals = tables[name]

            # Use the first row as reference.
            non_zeros = list(vals[0].nonzero()[0])

            # If all columns in the first row are non zero, there's no point
            # in continuing.
            if len(non_zeros) == shape(vals)[1]:
                continue

            # If we only have one row (IP) we just need the nonzero columns.
            if shape(vals)[0] == 1:
                if list(non_zeros):
                    non_zeros.sort()
                    non_zero_columns[name] = (i, non_zeros)

                    # Compress values.
                    tables[name] = vals[:, non_zeros]
                    i += 1

            # Check if the remaining rows are nonzero in the same positions, else expand.
            else:
                for j in range(shape(vals)[0] - 1):
                    # All rows must have the same non-zero columns
                    # for the optimization to work (at this stage).
                    new_non_zeros = list(vals[j+1].nonzero()[0])
                    if non_zeros != new_non_zeros:
                        non_zeros = non_zeros + [c for c in new_non_zeros if not c in non_zeros]
                        # If this results in all columns being non-zero, continue.
                        if len(non_zeros) == shape(vals)[1]:
                            continue

                # Only add nonzeros if it results in a reduction of columns.
                if len(non_zeros) != shape(vals)[1]:
                    if list(non_zeros):
                        non_zeros.sort()
                        non_zero_columns[name] = (i, non_zeros)

                        # Compress values.
                        tables[name] = vals[:, non_zeros]
                        i += 1

    # Check if we have some zeros in the tables.
    names_zeros = contains_zeros(tables, format_epsilon)

    # Add non-zero column, zero and ones info to inverse_name_map
    # (so we only need to pass around one name_map to code generating functions).
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

    # If we found non zero columns we might be able to reduce number of tables further.
    if non_zero_columns:
        # Try reducing the tables. This is possible if some tables have become
        # identical as a consequence of compressing the tables.
        # This happens with e.g., gradients of linear basis
        # FE0 = {-1,0,1}, nzc0 = [0,2]
        # FE1 = {-1,1,0}, nzc1 = [0,1]  -> FE0 = {-1,1}, nzc0 = [0,2], nzc1 = [0,1].

        # Call old utility function again.
        nm, inv_nm = unique_tables(tables, format_epsilon)

        # Update name maps.
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

        # Get new names of tables with all ones (for vector constants).
        names = get_ones(tables, format_epsilon)

        # Because these tables now contain ones as a consequence of compression
        # we still need to consider the non-zero columns when looking up values
        # in coefficient arrays. The psi entries can however we neglected and we
        # don't need to tabulate the values (if option is set).
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

    return (inverse_name_map, tables)

def unique_tables(tables, format_epsilon):
    """Removes tables with redundant values and returns a name_map and a
    inverse_name_map. E.g.,

    tables = {a:[0,1,2], b:[0,2,3], c:[0,1,2], d:[0,1,2]}
    results in:
    tables = {a:[0,1,2], b:[0,2,3]}
    name_map = {a:[c,d]}
    inverse_name_map = {a:a, b:b, c:a, d:a}."""

    name_map = {}
    inverse_name_map = {}
    names = sorted(tables.keys())
    mapped = []

    # Loop all tables to see if some are redundant.
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

            # Check if dimensions match.
            if shape(val0) == shape(val1):
                # Check if values are the same.
                if sqrt(((val0 - val1)*(val0 - val1)).sum()) < format_epsilon:
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
        if not name in inverse_name_map:
            inverse_name_map[name] = name

    return (name_map, inverse_name_map)

def get_ones(tables, format_epsilon):
    "Return names of tables for which all values are 1.0."
    names = []
    for name in tables:
        vals = tables[name]
        one = True
        for r in range(shape(vals)[0]):
            for c in range(shape(vals)[1]):
                if abs(vals[r][c] - 1.0) > format_epsilon:
                    one = False
        if one:
            names.append(name)
    return names

def contains_zeros(tables, format_epsilon):
    "Checks if any tables contains only zeros."

    names = []
    for name in tables:
        vals = tables[name]
        zero = True
        for r in range(shape(vals)[0]):
            # Loop as long as we have zeros.
            if not zero:
                break
            for c in range(shape(vals)[1]):
                # If just one value is different from zero, break loops.
                if abs(vals[r][c]) > format_epsilon:
                    # If we find a non-zero value break loop.
                    zero = False
                    break
        if zero:
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
                ffc_assert(tuple(key0 + key1) not in new, "This is not supposed to happen.")
                new[tuple(key0 + key1)] = val0 + val1

        return new

    # Create permutations by calling this function recursively.
    # This is only used for rank > 2 tensors I think.
    if len(expr) > 2:
        new = permutations(expr[0:2])
        return permutations(new + expr[2:])
