"Code generator for quadrature representation"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2007-03-16 -- 2007-06-01"
__copyright__ = "Copyright (C) 2007 Kristian B. Oelgaard"
__license__  = "GNU GPL Version 2"

# Modified by Anders Logg 2007

# Python modules
import os

# FFC common modules
#from ffc.common.constants import *
#from ffc.common.utils import *

# FFC language modules
from ffc.compiler.language.index import *
from ffc.compiler.language.restriction import *

# FFC code generation modules
from ffc.compiler.codegeneration.common.codegenerator import *
#from ffc.compiler.codegeneration.common.utils import *

# FFC tensor representation modules
#from ffc.compiler.representation.tensor.multiindex import *

# Should be in dictionary!!
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

def generate_psi_declaration(tensor_number, psi_indices, vindex, aindices, bindices,\
                                   num_quadrature_points, num_dofs, format):

    # Prefetch formats to speed up code generation
    format_secondary_index  = format["secondary index"]

    multi_index = []

    indices = ""
    for index in psi_indices:
        if index == vindex:
            indices = format_secondary_index(index_names[index.type](index.index)) + indices
        else:
            indices += format_secondary_index(index_names[index.type](index([], aindices, bindices, [])))
            multi_index += [index([], aindices, bindices, [])]

    name = format["table declaration"] + format["psis"] + format_secondary_index("t%d" %tensor_number)\
              + indices + format["matrix access"](num_quadrature_points, num_dofs)

    return (name, multi_index)

def generate_psi_entry(tensor_number, aindices, bindices, psi_indices, vindices, format):

    # Prefetch formats to speed up code generation
    format_secondary_index  = format["secondary index"]

    primary_indices = [format["first free index"], format["second free index"]]

    indices = ""
    for index in psi_indices:
        if index in vindices:
            indices = format_secondary_index(index_names[index.type](index.index)) + indices
            dof_num = index(primary_indices, aindices, bindices, [])
        else:
            indices += format_secondary_index(index_names[index.type](index([], aindices, bindices, [])))

    entry = format["psis"] + format_secondary_index("t%d" %tensor_number)\
              + indices + format["matrix access"](format["integration points"], dof_num)

    return entry

def generate_loop(name, value, loop_vars, Indent, format, connect = None):
    "This function generates a loop over a vector or matrix."

    code = []

    # Prefetch formats to speed up code generation
    format_loop     = format["loop"]

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

    group_tensors = {}
    for i in range(len(tensors)):
        tens = tensors[i]
        num_points = len(tens.quadrature.weights)
        if num_points in group_tensors:
            group_tensors[num_points] += [i]
        else:
            group_tensors[num_points] = [i]

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

def generate_signs(tensors, format):
    "Generate list of declarations for computation of signs"
    code = []
    for j in range(len(tensors)):
        monomial = tensors[j].monomial
        # Inspect each basis function (identified by its index)
        # and check whether sign changes are relevant.
        for basisfunction in monomial.basisfunctions:
            index = basisfunction.index
            values = []
            necessary = False
            element = basisfunction.element
            declarations = []
            dof_entities = DofMap(element).dof_entities();

            # Go through the topological entities associated
            # with each basis function/dof. If the element is
            # a piola mapped element and the basis function is
            # associated with an edge, we calculate the
            # possible sign change.
            for no in dof_entities:
                (entity, entity_no) = dof_entities[no]
                if entity == 1 and element.space_mapping(no) == Mapping.PIOLA:
                    necessary = True
                    values += [format["call edge sign"](entity_no)]
                else:
                    values += ["1"]
            num_dofs = str(len(dof_entities))
            name = format["sign tensor declaration"](format["signs"] + format["secondary index"]\
                   (str(j))+ format["secondary index"](str(index.index)) + format["array access"](num_dofs))
            value = format["block"](format["separator"].join([str(val) for val in values]))
 
            # Add declarations for this basis function to the code
            code += [(name,value)]
                    
    if necessary:
        code.insert(0, format["comment"]("Compute signs"))
        code.insert(0, format["snippet edge signs"](2))
        return (code, True)
    else:
        return ([], False) # Return [] is the case of no sign changes...)

def add_sign(value, j, i, format):
    if value:
        value = format["grouping"](value)
        for k in range(len(i)):
            value = format["multiply"]([format["signs"] + format["secondary index"]\
                   (str(j))+ format["secondary index"](str(k)) + format["array access"](i[k]),value])
    return value

def generate_factor(tensor, a, bgindices, format):
    "Optimise level 2 and 3"

# From tensorgenerator    
        # Compute product of factors outside sum
    factors = []
    for j in range(len(tensor.coefficients)):
        c = tensor.coefficients[j]
        if not c.index.type == Index.AUXILIARY_G:
            offset = tensor.coefficient_offsets[c]
            coefficient = format["coefficient"](c.n1.index, c.index([], a, [], [])+offset)
            for l in range(len(c.ops)):
                op = c.ops[len(c.ops) - 1 - l]
                if op == Operators.INVERSE:
                    coefficient = format["inverse"](coefficient)
                elif op == Operators.ABS:
                    coefficient = format["absolute value"](coefficient)
                elif op == Operators.SQRT:
                    coefficient = format["sqrt"](coefficient)
            factors += [coefficient]
    for t in tensor.transforms:
        if not (t.index0.type == Index.AUXILIARY_G or  t.index1.type == Index.AUXILIARY_G):
            factors += [format["transform"](t.type, t.index0([], a, [], []), \
                                                    t.index1([], a, [], []), \
                                                    t.restriction),]
    monomial = format["multiply"](factors)
    if monomial:
        f_out = [monomial]
    else:
        f_out = []
    
    # Compute sum of monomials inside sum
    terms = []
    for b in bgindices:
        factors = []
        for j in range(len(tensor.coefficients)):
            c = tensor.coefficients[j]
            if c.index.type == Index.AUXILIARY_G:
                offset = tensor.coefficient_offsets[c]
                coefficient = format["coefficient"](c.n1.index, c.index([], a, [], b)+offset)
                for l in range(len(c.ops)):
                    op = c.ops[len(c.ops) - 1 - l]
                    if op == Operators.INVERSE:
                        coefficient = format["inverse"](coefficient)
                    elif op == Operators.ABS:
                        coefficient = format["absolute value"](coefficient)
                    elif op == Operators.SQRT:
                        coefficient = format["sqrt"](coefficient)
                factors += [coefficient]
        for t in tensor.transforms:
            if t.index0.type == Index.AUXILIARY_G or t.index1.type == Index.AUXILIARY_G:
                factors += [format["transform"](t.type, t.index0([], a, [], b), \
                                                        t.index1([], a, [], b), \
                                                        t.restriction)]
        terms += [format["multiply"](factors)]

    f_in = format["add"](terms)
    if f_in: f_in = [format["grouping"](f_in)]
    else: f_in = []

#        f_tot = format["multiply"](f_out + f_in)

    return (f_out, f_in)

def generate_factor_old(tensor, a, b, format):
    "Optimise level 0 and 1"
    # Compute product of factors outside sum
    factors = []
    for j in range(len(tensor.coefficients)):
        c = tensor.coefficients[j]
        if not c.index.type == Index.AUXILIARY_G:
            offset = tensor.coefficient_offsets[c]
            coefficient = format["coefficient"](c.n1.index, c.index([], a, [], [])+offset)
            for l in range(len(c.ops)):
                op = c.ops[len(c.ops) - 1 - l]
                if op == Operators.INVERSE:
                    coefficient = format["inverse"](coefficient)
                elif op == Operators.ABS:
                    coefficient = format["absolute value"](coefficient)
                elif op == Operators.SQRT:
                    coefficient = format["sqrt"](coefficient)
            factors += [coefficient]
    for t in tensor.transforms:
        if not (t.index0.type == Index.AUXILIARY_G or  t.index1.type == Index.AUXILIARY_G):
            factors += [format["transform"](t.type, t.index0([], a, [], []), \
                                                    t.index1([], a, [], []), \
                                                    t.restriction),]
    # Compute sum of monomials inside sum
    for j in range(len(tensor.coefficients)):
        c = tensor.coefficients[j]
        if c.index.type == Index.AUXILIARY_G:
            offset = tensor.coefficient_offsets[c]
            coefficient = format["coefficient"](c.n1.index, c.index([], a, [], b)+offset)
            for l in range(len(c.ops)):
                op = c.ops[len(c.ops) - 1 - l]
                if op == Operators.INVERSE:
                    coefficient = format["inverse"](coefficient)
                elif op == Operators.ABS:
                    coefficient = format["absolute value"](coefficient)
                elif op == Operators.SQRT:
                    coefficient = format["sqrt"](coefficient)
            factors += [coefficient]
    for t in tensor.transforms:
        if t.index0.type == Index.AUXILIARY_G or t.index1.type == Index.AUXILIARY_G:
            factors += [format["transform"](t.type, t.index0([], a, [], b), \
                                                        t.index1([], a, [], b), \
                                                        t.restriction)]
    if tensor.determinant:
        d0 = format["power"](format["determinant"], tensor.determinant)
        d = format["multiply"]([format["scale factor"], d0])
    else:
        d = format["scale factor"]
    return [format["multiply"](factors + [d])]

def generate_factor_old2(tensor, a, b, format):
    "Optimise level 0 and 1"
    # Compute product of factors outside sum
    factors = []
    for j in range(len(tensor.coefficients)):
        c = tensor.coefficients[j]
        if not c.index.type == Index.AUXILIARY_G:
            offset = tensor.coefficient_offsets[c]
            if offset:
                coefficient = format["coeff"] + format["matrix access"](c.n1.index,\
                              format["add"]([c.index([], a, [], []), str(offset)]))
            else:
                coefficient = format["coeff"] + format["matrix access"](c.n1.index, c.index([], a, [], []))
            for l in range(len(c.ops)):
                op = c.ops[len(c.ops) - 1 - l]
                if op == Operators.INVERSE:
                    coefficient = format["inverse"](coefficient)
                elif op == Operators.ABS:
                    coefficient = format["absolute value"](coefficient)
                elif op == Operators.SQRT:
                    coefficient = format["sqrt"](coefficient)
            factors += [coefficient]
    for t in tensor.transforms:
        if not (t.index0.type == Index.AUXILIARY_G or  t.index1.type == Index.AUXILIARY_G):
            factors += [format["transform"](t.type, t.index0([], a, [], []), \
                                                    t.index1([], a, [], []), \
                                                    t.restriction),]
    # Compute sum of monomials inside sum
    for j in range(len(tensor.coefficients)):
        c = tensor.coefficients[j]
        if c.index.type == Index.AUXILIARY_G:
            offset = tensor.coefficient_offsets[c]
            if offset:
                coefficient = format["coeff"] + format["matrix access"](c.n1.index,\
                              format["add"]([c.index([], a, [], b), str(offset)]))
            else:
                coefficient = format["coeff"] + format["matrix access"](c.n1.index, c.index([], a, [], b))

            for l in range(len(c.ops)):
                op = c.ops[len(c.ops) - 1 - l]
                if op == Operators.INVERSE:
                    coefficient = format["inverse"](coefficient)
                elif op == Operators.ABS:
                    coefficient = format["absolute value"](coefficient)
                elif op == Operators.SQRT:
                    coefficient = format["sqrt"](coefficient)
            factors += [coefficient]
    for t in tensor.transforms:
        if t.index0.type == Index.AUXILIARY_G or t.index1.type == Index.AUXILIARY_G:
            factors += [format["transform"](t.type, t.index0([], a, [], b), \
                                                        t.index1([], a, [], b), \
                                                        t.restriction)]
    if tensor.determinant:
        d0 = format["power"](format["determinant"], tensor.determinant)
        d = format["multiply"]([format["scale factor"], d0])
    else:
        d = format["scale factor"]
    return [format["multiply"](factors + [d])]

def generate_factor3(tensor, a, bgindices, format):
    "Optimise level 3"

# From tensorgenerator    
        # Compute product of factors outside sum
    factors = []
    for j in range(len(tensor.coefficients)):
        c = tensor.coefficients[j]
        if not c.index.type == Index.AUXILIARY_G:
            offset = tensor.coefficient_offsets[c]
            if offset:
                coefficient = format["coeff"] + format["matrix access"](str(c.n1.index),\
                              format["add"]([str(c.index([], a, [], [])), str(offset)]))
            else:
                coefficient = format["coeff"] + format["matrix access"](c.n1.index, c.index([], a, [], []))
            for l in range(len(c.ops)):
                op = c.ops[len(c.ops) - 1 - l]
                if op == Operators.INVERSE:
                    coefficient = format["inverse"](coefficient)
                elif op == Operators.ABS:
                    coefficient = format["absolute value"](coefficient)
                elif op == Operators.SQRT:
                    coefficient = format["sqrt"](coefficient)
            factors += [coefficient]
    for t in tensor.transforms:
        if not (t.index0.type == Index.AUXILIARY_G or  t.index1.type == Index.AUXILIARY_G):
            factors += [format["transform"](t.type, t.index0([], a, [], []), \
                                                    t.index1([], a, [], []), \
                                                    t.restriction),]
    monomial = format["multiply"](factors)
    if monomial:
        f_out = [monomial]
    else:
        f_out = []
    
    # Compute sum of monomials inside sum
    terms = []
    for b in bgindices:
        factors = []
        for j in range(len(tensor.coefficients)):
            c = tensor.coefficients[j]
            if c.index.type == Index.AUXILIARY_G:
                offset = tensor.coefficient_offsets[c]
                if offset:
                    coefficient = format["coeff"] + format["matrix access"](str(c.n1.index),\
                                  format["add"]([str(c.index([], a, [], b)), str(offset)]))
                else:
                    coefficient = format["coeff"] + format["matrix access"](c.n1.index, c.index([], a, [], b))
                for l in range(len(c.ops)):
                    op = c.ops[len(c.ops) - 1 - l]
                    if op == Operators.INVERSE:
                        coefficient = format["inverse"](coefficient)
                    elif op == Operators.ABS:
                        coefficient = format["absolute value"](coefficient)
                    elif op == Operators.SQRT:
                        coefficient = format["sqrt"](coefficient)
                factors += [coefficient]
        for t in tensor.transforms:
            if t.index0.type == Index.AUXILIARY_G or t.index1.type == Index.AUXILIARY_G:
                factors += [format["transform"](t.type, t.index0([], a, [], b), \
                                                        t.index1([], a, [], b), \
                                                        t.restriction)]
        terms += [format["multiply"](factors)]

    f_in = format["add"](terms)
    if f_in: f_in = [format["grouping"](f_in)]
    else: f_in = []

#        f_tot = format["multiply"](f_out + f_in)

    return (f_out, f_in)


def values_level_0(indices, vindices, aindices, b0indices, bgindices, tensor, tensor_number, weight, format, name_map):
    # Generate value (expand multiplication - optimise level 0)
    format_multiply = format["multiply"]
    format_new_line = format["new line"]
    values = []
    for a in aindices:
        for b0 in b0indices:
            for bg in bgindices:
                factor = generate_factor_old(tensor, a, bg, format)
                values += [format_multiply([generate_psi_entry(tensor_number, a,\
                           b0, psi_indices, vindices, format) for psi_indices in indices] +\
                           weight + factor) + format_new_line]

    return (values, [])

def values_level_1(indices, vindices, aindices, b0indices, bgindices, tensor, tensor_number, weight, format, name_map):
    # Generate brackets of geometry and reference terms (terms outside sum are
    # multiplied with each term in the sum  - optimise level 1)
    format_multiply = format["multiply"]
    format_group = format["grouping"]
    format_add = format["add"]
    format_new_line = format["new line"]
    values = []
    for a in aindices:
        r = []
        for b0 in b0indices:
            r += [format_multiply([generate_psi_entry(tensor_number, a,\
                  b0, psi_indices, vindices, format) for psi_indices in indices] + weight)]

        if 1 < len(r):
            ref = format_group(format_add(r))
        else:
            ref = r[0]

        geo = [format_multiply(generate_factor_old(tensor, a, bg, format)) for bg in bgindices]
        if 1 < len(geo):
            geo = format_group(format_add(geo))
        else:
            geo = geo[0]
        values += [format_multiply([ref,geo]) + format_new_line]

    return (values, [])

def values_level_2(indices, vindices, aindices, b0indices, bgindices, tensor, tensor_number, weight, format, name_map):
    # Generate brackets of geometry and reference terms, distinguish between terms outside and 
    # inside sum (geometry only). - optimise level 2)

    format_multiply = format["multiply"]
    format_group = format["grouping"]
    format_add = format["add"]
    format_new_line = format["new line"]

    values = []

    for a in aindices:
        r = []
        for b0 in b0indices:
            r += [format_multiply([generate_psi_entry(tensor_number, a,\
                  b0, psi_indices, vindices, format) for psi_indices in indices] + weight)]

        if 1 < len(r):
            ref = format_group(format_add(r))
        else:
            ref = r[0]

        # Get geometry terms from inside sum, and outside sum
        geo_out, geo_in = generate_factor(tensor, a, bgindices, format)
        if 1 < len(geo_in):
            geo_in = [format_group(format_add(geo_in))]

        if tensor.determinant:
            d0 = format["power"](format["determinant"], tensor.determinant)
            d = [format["multiply"]([format["scale factor"], d0])]
        else:
            d = [format["scale factor"]]

        geo = format_multiply(geo_out + geo_in + d)
        values += [format_multiply([ref,geo]) + format_new_line]
    return (values, [])

def values_level_3(indices, vindices, aindices, b0indices, bgindices, tensor, tensor_number, weight, format, name_map):
    # Generate value - optimise level 3 (Based on level 2 but use loops to multiply reference and geo tensor)

    format_multiply = format["multiply"]
    format_group    = format["grouping"]
    format_add      = format["add"]
    format_new_line = format["new line"]
    format_space    = format["space"]
    format_psis     = format["psis"]

    # Get list of free secondary loop indices
    list_indices = format["free secondary indices"]
    values = []
    vals = []
    secondary_loop = []
    sec_indices = []
    for index in vindices:
        if index.type == Index.SECONDARY and len(index.range) > 1:
            sec_indices += [index]

            # Generate loop variables
            old_ind = [d for d in list_indices]
            m = len(old_ind)
            g = 0
            # If list of loop indices is not long enough generate some more
            while m - 1 < index.index:
                new_ind = [old_ind[i] + list_indices[j] for i in range(g, len(old_ind))\
                                                      for j in range(len(list_indices))]
                g = len(new_ind)
                old_ind += new_ind
                m = len(old_ind)
#            print "new_ind: ", new_ind
#            print "index: ", index.index
            # Pick index and generate information for loop generation
            index_name = old_ind[index.index]
            secondary_loop += [[index_name, 0, len(index.range)]]

#    print "sec_indices: ", sec_indices
    for a in aindices:
        # Change secondary index value to loop indices, for basis function indices
        for i in range(len(sec_indices)):
            a[sec_indices[i].index] = secondary_loop[i][0]
        r = []
        for b0 in b0indices:
            r += [format_multiply([generate_psi_entry(tensor_number, a,\
                  b0, psi_indices, vindices, format) for psi_indices in indices] + weight)]

        if 1 < len(r):
            ref = format_group(format_add(r))
        else:
            ref = r[0]

        # Get geometry terms from inside sum, and outside sum
        geo_out, geo_in = generate_factor3(tensor, a, bgindices, format)
        if 1 < len(geo_in):
            geo_in = [format_group(format_add(geo_in))]

        if tensor.determinant:
            d0 = format["power"](format["determinant"], tensor.determinant)
            d = [format["multiply"]([format["scale factor"], d0])]
        else:
            d = [format["scale factor"]]

        geo = format_multiply(geo_out + geo_in + d)

        vals += [format_multiply([ref,geo]) + format_new_line]

    # Only use values that are unique
    for val in vals:
        if not val in values:
            values += [val]
#    print "\nvalues[0]: ", values[0]

#    print "name map: ", name_map

    # Map values to correct table
    if name_map:
        map_values = []
        for val in values:
            for name in name_map:
                val = val.replace(name, name_map[name])
            map_values += [val]
        values = map_values

    return (values, secondary_loop)

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
                                      a, b, num_quadrature_points, num_dofs, format)

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

def generate_load_table(tensors):
    "Generate header to load psi tables"

    function = """
void load_table(double A[][%d], const char* name, int m)
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

    try:
        file = open("load_table.h", "r")
        code = file.read()
        code = code.split("\n")
        file.close()
    except IOError:
        code = False
    if not code:
        # Create code
        code = "#include <fstream>"
        for dof in group_num_dofs:
            code += function %(dof,dof)
    else:
        for dof in group_num_dofs:
            new_lines = function %(dof,dof)
            new_lines = new_lines.split("\n")
            if not new_lines[1] in code:
                code += new_lines

        code = "\n".join(code)
    # Write code
    file = open("load_table.h", "w")
    file.write(code)
    file.close()

def save_psis(tensors, facet0, facet1, Indent, format, tables):
    "Save psis tables instead of tabulating them"

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

    # Get list of psis
#    psis = tensor.Psis
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
#    outputs = []

    if tables:
#        print "save tables"

        for names in tables:
#            print "names: ", names
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
#            outputs += [out %(name, num_quadrature_points, num_dofs, name)]
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
#    code += [Indent.indent(load_int %("load") + format_end_line)]
    code += [Indent.indent(format["if"] + format["grouping"]("load"))]
    code += [Indent.indent(format["block begin"])]
    code += [Indent.indent("std::cout << \"Loading tables\" << std::endl;")]

    for load in loads:
        code += [Indent.indent(load + format_end_line)]

    code += [(Indent.indent("load"),"0")]
#    code += [Indent.indent("save();")]
    code += [Indent.indent(format["block end"])]
#    for output in outputs:
#        code += [Indent.indent(output)]

#    file = open("load.table", "w")
#    file.write("0")
#    file.close

    return code

def unique_tables(tensors, format):
    # Determine if some tensors have the same tables (and same names)

    name_map = {}
    tables = {}

    for tensor_number in range(len(tensors)):

        tensor = tensors[tensor_number]
        table = get_names_tables(tensor, tensor_number, format)
        for name in table:
            tables[name] = table[name]

#    print "tables 2: ", tables
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
#                    print "numpy.shape(val0): ", numpy.shape(val0)
#                    print "numpy.shape(val1): ", numpy.shape(val1)
#                    print "val0: ", val0
#                    print "val1: ", val1
                    diff = val1 - val0
                    # Check if values are the same
                    if not diff.any():
                        if name0 in name_map:
                            name_map[name0] += [name1]
                            if name1 in name_map:
                                del name_map[name1]

    inverse_name_map = {}
    for name in name_map:
#        print "name: ", name
        # Get name of psi
        name_strip = name.split()[-1].split("[")[0]
#        print "name_strip: ", name_strip

        maps = name_map[name]
        for m in maps:
#            print "m: ", m
            # Get name of psi
            m_strip = m.split()[-1].split("[")[0]
#            print "m_strip: ", m_strip
            inverse_name_map[m_strip] = name_strip
            del tables[m]
#    print "inverse_name_map: ", inverse_name_map

    return (inverse_name_map, tables)



