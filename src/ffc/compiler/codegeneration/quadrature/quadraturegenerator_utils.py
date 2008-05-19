"Utility functions for quadrature representation"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2007-03-16 -- 2008-02-05"
__copyright__ = "Copyright (C) 2007-2008 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Anders Logg 2007

# Python modules
import os
from sets import Set

# FFC language modules
from ffc.compiler.language.index import *
from ffc.compiler.language.restriction import *

# FFC code generation modules
from ffc.compiler.codegeneration.common.codegenerator import *

# FIXME: This should be in dictionary!!
index_names = {0: lambda i: "f%s" %(i), 1: lambda i: "p%s" %(i), 2: lambda i: "s%s" %(i),\
               4: lambda i: "fu%s" %(i), 5: lambda i: "pj%s" %(i), 6: lambda i: "c%s" %(i), 7: lambda i: "a%s" %(i)}

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

def add_to_dict(dictionary, keys, val, add_if_value_exists = True):

    # If this is the last key add value to dictionary
    if len(keys) == 1:
        if keys[0] in dictionary:
            values = dictionary[keys[0]]
            if val[0] in values and add_if_value_exists:
                values += val
            elif val[0] not in values:
                values += val
            dictionary[keys[0]] = values
        else:
            dictionary[keys[0]] = val
    else:
        if not keys[0] in dictionary:
            dictionary[keys[0]] = {}
        add_to_dict(dictionary[keys[0]], keys[1:], val, add_if_value_exists)

def get_dict_keys_vals(v, keys):

    keys_vals = []
    for key in v:
        val = v[key]
        if isinstance(val, dict):
            keys_vals += get_dict_keys_vals(val, keys + [key])
        else:
            k = keys + [key]
            keys_vals.append((k, val))
    return keys_vals


def generate_psi_declaration(tensor_number, psi_indices, vindex, aindices, bindices,\
                                   num_quadrature_points, num_dofs, format, qeindices = []):

    # Prefetch formats to speed up code generation
    format_secondary_index  = format["secondary index"]
    format_table            = format["table declaration"]
    format_psis             = format["psis"]
    format_matrix_access    = format["matrix access"]
    multi_index = []

    indices = ""
    for index in psi_indices:
        # We only need QE psi tables if they're primary indices
        if index in qeindices and not index.type == Index.PRIMARY:
            return (None, None)
        if index == vindex:
            indices = format_secondary_index(index_names[index.type](index.index)) + indices
        else:
            indices += format_secondary_index(index_names[index.type](index([], aindices, bindices, [])))
            multi_index += [index([], aindices, bindices, [])]

    name = format_table + format_psis + format_secondary_index("t%d" %tensor_number)\
              + indices + format_matrix_access(num_quadrature_points, num_dofs)

    return (name, multi_index)

def generate_psi_entry(tensor_number, aindices, bindices, psi_indices, vindices, name_map, format, qeindices = []):

    # Prefetch formats to speed up code generation
    primary_indices         = [format["first free index"], format["second free index"]]
    format_secondary_index  = format["secondary index"]
    format_psis             = format["psis"]
    format_matrix_access    = format["matrix access"]
    format_ip               = format["integration points"]

    indices = ""
    for index in psi_indices:
        if index in qeindices:
            return ""
        if index in vindices:
            indices = format_secondary_index(index_names[index.type](index.index)) + indices
            dof_num = index(primary_indices, aindices, bindices, [])
        else:
            indices += format_secondary_index(index_names[index.type](index([], aindices, bindices, [])))

    name = format_psis + format_secondary_index("t%d" %tensor_number) + indices
    if name in name_map:
        entry = name_map[name] + format_matrix_access(format_ip, dof_num)
    else:
        entry = name + format_matrix_access(format_ip, dof_num)

    return entry

def generate_psi_entry2(tensor_number, aindices, bindices, psi_indices, vindices, name_map, non_zero_columns, format, qeindices = []):

    # Prefetch formats to speed up code generation
    primary_indices         = [format["first free index"], format["second free index"]]
    format_secondary_index  = format["secondary index"]
    format_psis             = format["psis"]
    format_matrix_access    = format["matrix access"]
    format_ip               = format["integration points"]

    indices = ""
    dof_range = -1
    for index in psi_indices:
#        if index in qeindices:
#            return (["P", "F"], dof_range, "")
        if index in vindices:
            indices = format_secondary_index(index_names[index.type](index.index)) + indices
            loop_dof = index(primary_indices, aindices, bindices, [])
        else:
            indices += format_secondary_index(index_names[index.type](index([], aindices, bindices, [])))

    name = format_psis + format_secondary_index("t%d" %tensor_number) + indices

#    print "name: ", name
    # For QE, we only need a psi entry if we don't have a primary index
    for index in psi_indices:
#        print "index: ", index
#        print "qeindices: ", qeindices
        if index in qeindices and not index.type == Index.PRIMARY:
            name = ""
    if name:
        if name in name_map:
            name = name_map[name]
        dof_num = loop_dof
        if name in non_zero_columns:
            i, cols = non_zero_columns[name]
            dof_range = len(cols)
            if dof_range == 1:
                dof_num = format["nonzero columns"](i) + format["array access"]("0")
            else:
                dof_num = format["nonzero columns"](i) + format["array access"](loop_dof)
            entry = name + format_matrix_access(format_ip, dof_num)
        else:
            entry = name + format_matrix_access(format_ip, loop_dof)
            if not name:
                entry = name
    else:
        dof_num = loop_dof
        entry = name

    return ([loop_dof, dof_num], dof_range, entry)

def generate_loop(name, value, loop_vars, Indent, format, connect = None):
    "This function generates a loop over a vector or matrix."

    # Prefetch formats to speed up code generation
    format_loop             = format["loop"]
    code = []
    for ls in loop_vars:

        # Get index and lower and upper bounds
        index, lower, upper = (ls[0], ls[1], ls[2])

        # Loop index
        code += [Indent.indent(format_loop(index, lower, upper))]

        # Increase indentation
        Indent.increase()
        if name:
            # If this is the last loop, write values
            if index == loop_vars[-1][0]:
                if connect:
                    code += [connect(Indent.indent(name), value)]
                else:
                    code += [(Indent.indent(name), value)]

    # Decrease indentation
    if name:
        for i in range(len(loop_vars)):
            Indent.decrease()

    return code

def generate_loop2(lines, loop_vars, Indent, format):
    "This function generates a loop over a vector or matrix."

    # Prefetch formats to speed up code generation
    format_loop             = format["loop"]
    code = []
    for ls in loop_vars:
        # Get index and lower and upper bounds
        index, lower, upper = (ls[0], ls[1], ls[2])

        # Loop index
        code += [Indent.indent(format_loop(index, lower, upper))]

        # Increase indentation
        Indent.increase()
        # If this is the last loop, write values
        if index == loop_vars[-1][0]:
            # Decrease indentation
            Indent.decrease()
            code += [Indent.indent(format["block begin"])]
            # Increase indentation
            Indent.increase()
            code += [Indent.indent(l) for l in lines]
            # Decrease indentation
            Indent.decrease()
            code += [Indent.indent(format["block end"])]

    for i in range(len(loop_vars) - 1):
        Indent.decrease()

    return code

def generate_name_entry(loop_vars, format):

#    print loop_vars
    prim = loop_vars[0]
#    print "prim[0]: ", prim[0]
#    print "sec[0]: ", sec[0]
#    print "prim: ", prim
#    print "sec: ", sec
    # Test validity
    if not len(prim[0]) == len(prim[1]):
        raise RuntimeError, "Something is very wrong!!"

    # Generate name for entry in element tensor (move this to term generation??)
    name = ""
    rank = len(prim[0])
    if (rank == 0):
        # Entry is zero because functional is a scalar value
        entry = "0"
        # Generate name
        name =  format["element tensor quad"] + format["array access"](entry)
    elif (rank == 1):
        # Generate entry
        entry = prim[1][0]
        # Generate name
        name =  format["element tensor quad"] + format["array access"](entry)

    elif (rank == 2):
        entry = list(loop_vars[1])
        # Generate entry
        entry[0] = format["multiply"]([entry[0], str(prim[0][1])])
        name =  format["element tensor quad"] + format["array access"](format["add"](entry))
    else:
        raise RuntimeError, "Quadrature only support Functionals and Linear and Bilinear forms"
    return name

def generate_loop3(lines, loop_vars, Indent, format):
    "This function generates a loop over a vector or matrix."

    # Get primary and secondary loop info
#    print "loop_vars: ", loop_vars
    range_indices = loop_vars[0]
    indices       = loop_vars[1]
#    print "prim[0]: ", prim[0]
#    print "sec[0]: ", sec[0]
#    print "prim: ", prim
#    print "sec: ", sec
    # Test validity
    if not len(range_indices) == len(indices):
        raise RuntimeError, "Something is very wrong!!"

    new_loops = []
    for i in range(len(indices)):
        index = indices[i]
        lower = 0
        upper = range_indices[i]
        new_loops += [(index, lower, upper)]
#    print new_loops

    if new_loops:
        return generate_loop2(lines, new_loops, Indent, format)
    else:
        return lines


def extract_unique(aa, bb):
    "Remove redundant names and entries in the psi table"

    uaa = []
    ubb = []

    for i in range(len(aa)):
        a = aa[i]
        if not a in uaa:
            uaa += [a]
            ubb += [bb[i]]
    return (uaa, ubb)

def equal_loops(tensors):
    "Group tensor with an equal number of quadrature points and primary indices"

    # Loop all tensors and group number of quadrature points
    group_tensors = {}
    for i in range(len(tensors)):
        tens = tensors[i]
        num_points = len(tens.quadrature.weights)
        if num_points in group_tensors:
            group_tensors[num_points] += [i]
        else:
            group_tensors[num_points] = [i]

    # Then loop all groups of quadrature points and group after primary indices
    for g in group_tensors:
        tensor_nums = group_tensors[g]
        prims = {}
        for t in tensor_nums:
            tens = tensors[t]
            idims = tens.i.dims
            if tuple(idims) in prims:
                prims[tuple(idims)] += [t]
            else:
                prims[tuple(idims)] = [t]
        group_tensors[g] = prims

    return group_tensors

def get_names_tables(tensor, tensor_number, format):
    "Tabulate values of basis functions and their derivatives at quadrature points"

    tables = {}

    # Get list of psis
    psis = tensor.Psis

    # Loop psis
    for psi_number in range(len(psis)):

        # Get psi
        psi = psis[psi_number]

        # Get values of psi, list of indices and index of basisfunction for psi
        values, indices, vindex = psi[0], psi[1], psi[2]

        # Get the number of dofs
        num_dofs = len(vindex.range)

        # Get number of quadrature points
        num_quadrature_points = len(tensor.quadrature.weights)

        # Get lists of secondary indices and auxiliary_0 indices
        aindices = tensor.a.indices
        b0indices = tensor.b0.indices

        names = []
        multi_indices = []

        # Loop secondary and auxiliary indices to generate names and entries in the psi table
        for a in aindices:
            for b in b0indices:

                (name, multi_index) = generate_psi_declaration(tensor_number, indices, vindex,\
                                      a, b, num_quadrature_points, num_dofs, format, tensor.qei)
#                (name, multi_index) = generate_psi_declaration(tensor_number, indices, vindex,\
#                                      a, b, num_quadrature_points, num_dofs, format)

                if name:
                    names += [name]
                    multi_indices += [multi_index]

        # Remove redundant names and entries in the psi table
        names, multi_indices = extract_unique(names, multi_indices)

        # Loop names and tabulate psis
        for i in range(len(names)):

            # Get name
            name = names[i]

            # Get values from psi tensor, should have format values[dofs][quad_points]
            vals = values[tuple(multi_indices[i])]

            # Check if the values have the correct dimensions, otherwise quadrature is not correctly
            # implemented for the given form!!
            if numpy.shape(vals) != (num_dofs, num_quadrature_points):
                raise RuntimeError, "Quadrature is not correctly implemented for the given form!"

            # Generate values (FIAT returns [dof, quad_points] transpose to [quad_points, dof])
            value = numpy.transpose(vals)
            tables[name] = value

    return tables

def unique_psi_tables(tensors, optimisation_level, format):
    "Determine if some tensors have the same tables (and same names)"

    name_map = {}
    tables = {}
    non_zero_columns = {}

    # Loop tensors and get all tables
    for tensor_number in range(len(tensors)):

        tensor = tensors[tensor_number]
        table = get_names_tables(tensor, tensor_number, format)
        # Copy table to global dictionary of tables
        for name in table:
            tables[name] = table[name]

    # Initialise name map
    for name in tables:
        name_map[name] = []

    # Loop all tables to see if some are redundant
    for name0 in tables:
        val0 = numpy.array(tables[name0])
        for name1 in tables:
            # Don't compare values with self
            if not name0 == name1:
                val1 = numpy.array(tables[name1])
                # Check if dimensions match
                if numpy.shape(val0) == numpy.shape(val1):
                    # Compute difference
                    diff = val1 - val0
                    # Check if values are the same
                    if not diff.any():
                        if name0 in name_map:
                            name_map[name0] += [name1]
                            if name1 in name_map:
                                del name_map[name1]

    # Construc the inverse of the name map
    inverse_name_map = {}
    for name in name_map:
        # Get name of psi
        # FIXME: this is C++ dependent
        name_strip = name.split()[-1].split("[")[0]

        maps = name_map[name]
        for m in maps:
            # Get name of psi
            # FIXME: this is C++ dependent
            m_strip = m.split()[-1].split("[")[0]
            inverse_name_map[m_strip] = name_strip
            del tables[m]

    # If we have all ones we don't need this table
    if optimisation_level >= 10:
        ones = []
        for name in tables:
            vals = tables[name]
            one = True
            for r in range(numpy.shape(vals)[0]):
                for c in range(numpy.shape(vals)[1]):
                    if not vals[r][c] == 1.0:
                        one = False
            if one:
                ones.append(name)
        for one in ones:
              o_strip = one.split()[-1].split("[")[0]
              del tables[one]
#              print "one: "
#              print one
#              print name_map
              maps = name_map[one]
#              print "maps: ", maps
              inverse_name_map[o_strip] = ""
              for m in maps:
                  m_strip = m.split()[-1].split("[")[0]
                  inverse_name_map[m_strip] = ""
#              if not maps:
#                  inverse_name_map[o_strip] = ""
#              else:
#                  for m in maps:
#                      m_strip = m.split()[-1].split("[")[0]
#                      inverse_name_map[m_strip] = ""

#    print "inverse name map: ", inverse_name_map

    if optimisation_level <= 5:
        return (inverse_name_map, tables, non_zero_columns)
    else:
        # Extract the column numbers that are non-zero
        i = 0
        for name in tables:
            # Get values and take the first row as reference
            vals = tables[name]
#            print "vals: ", vals
            # Set values to zero if they are lower than threshold
            for r in range(numpy.shape(vals)[0]):
                for c in range(numpy.shape(vals)[1]):
                    if abs(vals[r][c]) < format["epsilon"]:
                        vals[r][c] = 0
#            print "vals: ", vals

#            print numpy.shape(vals)
#            print "vals[0].nozero(): ", vals[0].nonzero()[0]
            non_zeros = list(vals[0].nonzero()[0])
            # If all columns in the first row are non zero, there's no point
            # in continuing
            if len(non_zeros) == numpy.shape(vals)[1]:
                # Values in all positions, do not create map
#                non_zeros = range(numpy.shape(vals)[1])
#                non_zero_columns[name] = (i, non_zeros)
#                i += 1
                continue
#            print "vals[1]: ", vals[1]
#            print "vals[1].nozero(): ", vals[1].nonzero()

            # If we only have one row (IP) we just need the nonzero columns
            if numpy.shape(vals)[0] == 1:
                if list(non_zeros):
                    non_zeros.sort()
                    non_zero_columns[name] = (i, non_zeros)
                    i += 1
            # Check if the remaining rows are nonzero in the same positions
            else:
#                print "more than one row"
                for j in range(numpy.shape(vals)[0] - 1):
#                    print "row: ", j
                    # All rows must have the same non-zero columns
                    # for the optimization to work (at this stage)
                    new_non_zeros = list(vals[j+1].nonzero()[0])
                    if non_zeros != new_non_zeros:
                        non_zeros = non_zeros + [c for c in new_non_zeros if not c in non_zeros]
#                        print "\n\nno non zeros: ", name
#                        print "j: ", j+1
#                        print non_zeros
#                        print list(vals[j+1].nonzero()[0])
                        # If not all rows have the same non-zero columns, return
                        # all columns (assume no non-zeros)
#                        non_zeros = range(numpy.shape(vals)[1])
#                        print vals
#                        break
#                # Only add nonzeros if all rows were identical
                if list(non_zeros):
                    non_zeros.sort()
                    non_zero_columns[name] = (i, non_zeros)
                    i += 1
        return (inverse_name_map, tables, non_zero_columns)

def unique_weight_tables(tensors, format):
    "Determine if some tensors have the same tables (and same names)"

    name_map = {}
    tables = {}

    # Loop tensors and get all tables
    for tensor_number in range(len(tensors)):
        weights = tensors[tensor_number].quadrature.weights
#        name = tensor_number
#        name = format["table declaration"] + format["weights"](tensor_number, str(len(weights)))
        tables[tensor_number] = weights

#    print tables
    # Initialise name map
    for name in tables:
        name_map[name] = []

    # Loop all tables to see if some are redundant
    for name0 in tables:
        val0 = numpy.array(tables[name0])
        for name1 in tables:
            # Don't compare values with self
            if not name0 == name1:
                val1 = numpy.array(tables[name1])
#                print "val0: ", val0
#                print "val1: ", val1
                # Check if dimensions match
                if numpy.shape(val0) == numpy.shape(val1):
                    # Compute difference
                    diff = val1 - val0
#                    print "diff: ", diff
#                    print "diff.any(): ", diff.any()
                    # Check if values are the same
                    if not diff.any():
                        if name0 in name_map:
                            name_map[name0] += [name1]
                            if name1 in name_map:
                                del name_map[name1]
    return name_map

def generate_load_table(tensors):
    "Generate header to load psi tables"

    # FIXME: This function is VERY C++ and linux dependent
    function = """
void load_table(double A[][%d], const char* name, int m) const
{
  std::ifstream in(name, std::ios::in);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < %d; j++)
      in >> A[i][j];
  in.close();
}\n"""

#    load_int = """
#void load_table(unsigned int& load)
#{
#  std::ifstream in("load.table", std::ios::in);
#  in >> load;
#  in.close();
#}"""
#    save_int = """
#void save()
#{
#  std::ofstream out("load_f%d%d.table", std::ios::out);
#  out << 0;
#  out.close();
#}\n""" %(facet0, facet1)

    group_num_dofs = []

    for tensor in tensors:

        # Get list of psis
        psis = tensor.Psis

        # Loop psis
        for psi_number in range(len(psis)):

            # Get psi
            psi = psis[psi_number]

            # Get the number of dofs
            num_dofs = len(psi[2].range)

            if not num_dofs in group_num_dofs:
                group_num_dofs += [num_dofs]

    code = []
    for dof in group_num_dofs:
        new_lines = function %(dof,dof)
        new_lines = new_lines.split("\n")
        if not new_lines[1] in code:
            code += new_lines
    return "\n".join(code)

def save_psis(tensors, facet0, facet1, Indent, format, tables):
    "Save psis tables instead of tabulating them"

    # FIXME: This function is VERY C++ and linux dependent
    load_table = "load_table(%s, \"%s\", %d)"
    load_int = "load_table(%s)"
    out = """
std::cout.precision(17);
std::cout << "%s" << std::endl;
for (int j = 0; j < %d; j++)
{
  for (int k = 0; k < %d; k++)
    std::cout << %s[j][k] << " ";
  std::cout << std::endl;
}
std::cout << std::endl;
"""
    code = []
    format_floating_point = format["floating point"]
    format_epsilon        = format["epsilon"]
    format_end_line       = format["end line"]
    format_table          = format["table declaration"]
    format_float          = format["static float declaration"]
#    format_float          = format["const float declaration"]

    # Get list of psis
    # FIXME: C++ naming
    irank = "_i%d" %tensors[0].i.rank
    if (facet0 == None) and (facet1 == None):
        facets = ""
    elif (facet0 != None) and (facet1 == None):
        facets = "_f%d" %(facet0)
    else:
        facets = "_f%d%d" %(facet0, facet1)

    try:
        os.system("mkdir -p tables")

    except IOError:
        dirs = True

    loads = []

    if tables:
        for names in tables:
            # Get value and save as an array
            value = tables[names]
            shape = numpy.shape(value)
            name = names.split(" ")[-1].split("[")[0]
            a = []
            for j in range(shape[0]):
                for k in range(shape[1]):
                    val = value[j][k]
                    if abs(val) < format_epsilon:
                        val = 0.0
                    a += [format_floating_point(val)]
            a = " ".join(a)
            # FIXME: linux dependent
            file = open("tables/" + name + facets + irank + ".table", "w")
            file.write(a)
            file.close

            # FIXME: linux dependent
            # Generate code to load table
            loads += [load_table %(name, "tables/" + name + facets + irank + ".table", shape[0])]
            code += [Indent.indent(names.replace(format_table,format_float) + format_end_line)]
    else:
        for tensor_number in range(len(tensors)):
            tensor = tensors[tensor_number]
            tables = get_names_tables(tensor, tensor_number, format)

            for names in tables:
                # Get value and save as an array
                value = tables[names]
                shape = numpy.shape(value)
                name = names.split(" ")[-1].split("[")[0]
                a = []
                for j in range(shape[0]):
                    for k in range(shape[1]):
                        val = value[j][k]
                        if abs(val) < format_epsilon:
                            val = 0.0
                        a += [format_floating_point(val)]
                a = " ".join(a)
                file = open("tables/" + name + facets + irank + ".table", "w")
                file.write(a)
                file.close

                # Generate code to load table
                loads += [load_table %(name, "tables/" + name + facets + irank + ".table", shape[0])]
#                outputs += [out %(name, num_quadrature_points, num_dofs, name)]
                code += [Indent.indent(names.replace(format_table,format_float) + format_end_line)]

    code += [(Indent.indent(format["static uint declaration"] + "load"),"1")]
    code += [Indent.indent(format["if"] + format["grouping"]("load"))]
    code += [Indent.indent(format["block begin"])]

    load_code = [Indent.indent("std::cout << \"Loading tables\" << std::endl;")]
    for load in loads:
        code += [Indent.indent(load + format_end_line)]

    code += [(Indent.indent("load"),"0")]
    code += [Indent.indent(format["block end"])]

    return code


