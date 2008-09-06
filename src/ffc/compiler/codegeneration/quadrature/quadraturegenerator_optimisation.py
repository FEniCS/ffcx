"Optimisation functions for quadrature representation"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2008-02-05 -- 2008-02-05"
__copyright__ = "Copyright (C) 2008 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Python modules
import os
from sets import Set


# FFC language modules
from ffc.compiler.language.index import *
from ffc.compiler.language.restriction import *

# FFC code generation modules
from ffc.compiler.codegeneration.common.codegenerator import *

# Utility functions for quadraturegenerator
from quadraturegenerator_utils import *
import reduce_operations

def generate_code(raw_terms, geo_terms, optimise_level, Indent, format):

    format_comment      = format["comment"]
    format_add          = format["add"]
    format_add_equal    = format["add equal"]
    format_tensor       = format["element tensor quad"]
    format_array_access = format["array access"]
    format_G            = format["geometry tensor"]
    format_F            = format["function value"]
    format_float        = format["floating point"]
    format_float_decl   = format["float declaration"]

    exp_ops       = reduce_operations.expand_operations
    get_geo_terms = reduce_operations.get_geo_terms
    red_ops       = reduce_operations.reduce_operations
    count_ops     = reduce_operations.operation_count
    get_vars      = reduce_operations.get_variables

    terms_code = []
    total_num_ops = 0
    ip_terms = {}

    # A new dictionary for primary loop and entries
    terms = {}

    # Dictionary of functions that should be computed first for each
    # secondary loop before starting loops over primary indices.
    secondary_loops = {}
    functions = {}

    # Get loops and entries dictionary
    for loops, entry_dict in raw_terms.items():
#            print "\nLoops: ", loops
        p_loops, s_loops = loops
        # Create primary and secondary loops
        prim_loops = tuple([(p_loops[1][i], 0, p_loops[0][i]) for i in range(len(p_loops[0]))])
#            print "  prim_loops: ", prim_loops
        sec_loops = [(s_loops[1][i], 0, s_loops[0][i]) for i in range(len(s_loops[0]))]
#            print "  sec_loops: ", sec_loops

        # Loop entries and generate values
        for e, v in entry_dict.items():
#                print "\n    entry: ", e
            val = format_add(v)

            # Extract any constant terms, including geometry terms
            val = exp_ops(val, format)
            val = get_geo_terms(val, geo_terms, 0, format)

#                print "    val: ", val
            # Get any functions that we might have
            for sec in sec_loops:
#                    print "      Generating functions for sec_loop: ", sec
                val, used_vars = get_vars(val, functions, format, [sec[0]])
#                    print "\nUsed vars: ", used_vars
                if not sec in secondary_loops:
                    secondary_loops[sec] = Set(used_vars)
                else:
                    secondary_loops[sec] = secondary_loops[sec] | Set(used_vars)
#                    print "val: ", val

            # Add code to primary dictionary
            if prim_loops in terms:
                if e in terms[prim_loops]:
                    terms[prim_loops][e].append(val)
                else:
                    terms[prim_loops][e] = [val]
            else:
                terms[prim_loops] = {}
                terms[prim_loops][e] = [val]

#        print "functions: ", functions
#        print "secondary_loops: ", secondary_loops
    
    # Create list to declare function values
    function_declarations = []
    compute_functions = []
    if functions:
        function_declarations = ["", format_comment("Declare function values.")]
        function_declarations += [(format_float_decl + format_F + str(i),\
                            format_float(0)) for i in range(len(functions))]

        # Create inverse of functions dictionary
        inv_functions = dict([(v,k) for k,v in functions.items()])

        # Generate code to compute functions
        compute_functions = ["", format_comment("Compute function values.")]

        for sec, funcs in secondary_loops.items():
#                print "Funcs: ", funcs
            lines = []
            func_ops = 0
            for f in funcs:
                val = red_ops(inv_functions[f], format)
                func_ops += count_ops(val, format)
                lines.append(format_add_equal(f, val))

            # Compute count and add to global count
            func_ops *= sec[2]
            total_num_ops += func_ops
            compute_functions.append(format_comment("Number of operations to compute values = %d" % func_ops))
            compute_functions += generate_loop(lines, [sec], Indent, format)


    # Generate code for primary indices
    prim_lines = ["", format_comment("Loop primary indices.")]

    # Get primary loops and entries dictionary
    for prim_loops, entries in terms.items():

        # Code lines for all entries
        entry_lines = []
        entry_ops = 0
        # Get entries and values
        for e, v in entries.items():

            # Create value and count number of operations
            val = format_add(v)
            print "val: ", val

            num_ops = count_ops(val, format)
            entry_ops += num_ops

            # Generate comment, element tensor entry and add to lines
            entry_lines.append( format_comment\
                ("Number of operations to compute entry = %d" %(num_ops)) )
            entry_lines.append( format_add_equal( format_tensor + format_array_access(e) , val) )

        # Multiply number of operations by range of primary loops
        for p in prim_loops:
            entry_ops *= p[2]
        total_num_ops += entry_ops

        prim_lines.append(format_comment("Number of operations for primary indices = %d" % entry_ops))
        # Add code to primary loop code
        prim_lines += generate_loop(entry_lines, prim_loops, Indent, format)
        prim_lines.append("")

    # Add blocks of code
    terms_code += function_declarations + compute_functions + prim_lines

    return (terms_code, total_num_ops, ip_terms)

class Term:

    def __init__(self, tensor, list_indices):

        self.indices = [psi[1] for psi in tensor.Psis]
        self.vindices = [psi[2] for psi in tensor.Psis]

        # Get secondary indices
        self.secondary_indices = []
        for index in self.vindices:
            if index.type == Index.SECONDARY and len(index.range) > 1 and not index in tensor.qei:
                self.secondary_indices.append(index)

        # Generate loop variables
        sec_indices = self.secondary_indices  # Get list of secondary indices
        old_ind = [d for d in list_indices]   # Get copy of free secondary indices
        g = 0
        # If list of loop indices is not long enough generate some more
        while len(old_ind) - 1 < len(sec_indices):
            new_ind = [old_ind[i] + list_indices[j] for i in range(g, len(old_ind))\
                                                          for j in range(len(list_indices))]
            g = len(new_ind)
            old_ind += new_ind

        # Generate information for loop generation
        self.secondary_loops = [[old_ind[i], 0, len(sec_indices[i].range)] for i in range(len(sec_indices))]
#        print "secondary loops: ", self.secondary_loops

        irank = tensor.i.rank
        self.prim_dims = [0,]*irank
        self.prim_vars = [0,]*irank
        self.prim_dofs = [0,]*irank

        # Reduce aindices according to secondary loops
        aindices = tensor.a.indices
#        print "aindices: ", aindices
#        print "b0indices: ", b0indices
        # Reduce multi indices a and b0 (number of pertubations) according to current psi
        a_reduced = [[0],]*len(aindices[0])
#        print "a_reduced: ", a_reduced
        indices = []
        for i in self.indices: indices += i
        for index in indices:
            if index.type == Index.SECONDARY and not index in self.vindices:
                a_reduced[index.index] = range(len(index.range))
#                a_reduced[index.index] = index.range
        a_reduced = MultiIndex(a_reduced).indices
#        print "a_reduced: ", a_reduced
        self.reduced_aindices = a_reduced

def generate_factor(tensor, a, mapped_a, bgindices, b0, format, opt_level):
    "Generate factor"

    if not mapped_a:
        mapped_a = [index for index in a]

    format_multiply = format["multiply"]
    format_add      = format["add"]
    format_ip       = format["integration points"]
        # Compute product of factors outside sum
    trans_set = Set()
    factors = []
    factors_qe = []
    # If we have quadrature indices construct component access
    qe_access = {}
    if tensor.qei:
        # Substitute aindices for quadrature elements with direct component access
        for v in tensor.monomial.basisfunctions:
            if v.index in tensor.qei:

                # The number of components are equal to the number of subelements (for QE)
                num_comp = v.element.num_sub_elements()
                for i in range(num_comp):
                    elem = v.element.sub_element(i)
#                    print i
                    dim = elem.space_dimension()
                    if not i:
                        qe_access[i] = format_ip
                    else:
                        if dim > 1:
                            qe_access[i] = format_add([ "%d" %(dim*i), format_ip])
                        elif dim == 1:
                            qe_access[i] = format_add([ "%d" %i, format_ip])
                        else:
                            raise RuntimeError("Unexpected space_dimension!!")


#    print "qe_access: ", qe_access
#    print "b0indices: ", b0

    for j in range(len(tensor.coefficients)):
        c = tensor.coefficients[j]
#        print "c: ", c
#        print "c.index: ", c.index
#        print "c.n0.index: ", c.n0.index
#        print "c.n1.index: ", c.n1.index
#        print "c.e0: ", c.e0
#        print "c.e1: ", c.e1

        if not c.index.type == Index.AUXILIARY_G:
            offset = tensor.coefficient_offsets[c]
            access = str(c.index([], mapped_a, [], []))
            if c.index in tensor.qei:
                for v in tensor.monomial.basisfunctions:
                    if v.index == c.index:
#                        print "v.component: ", v.component
                        if v.component:
                            if len(v.component) == 1:
                                access = qe_access[v.component[0]([], mapped_a, b0, [])]
                            else:
                                raise RuntimeError("Error, more than one component index!!")
                        else:
                            access = qe_access[0]
#                            raise RuntimeError("Error!!")

            if offset:
                acs = format["add"]([access, str(offset)])
                try: acs = "%d" % eval(acs)
                except: pass
                coefficient = format["coeff"] + format["matrix access"](str(c.n1.index), acs)
            else:
                coefficient = format["coeff"] + format["matrix access"](c.n1.index, access)

            for l in range(len(c.ops)):
                op = c.ops[len(c.ops) - 1 - l]
                if op == Operators.INVERSE:
                    coefficient = format["inverse"](coefficient)
                elif op == Operators.MODULUS:
                    coefficient = format["absolute value"](coefficient)
                elif op == Operators.SQRT:
                    coefficient = format["sqrt"](coefficient)
            factors += [coefficient]
    for t in tensor.transforms:
        if not (t.index0.type == Index.AUXILIARY_G or  t.index1.type == Index.AUXILIARY_G):
            trans = format["transform"](t.type, t.index0([], a, [], []), \
                                                t.index1([], a, [], []), t.restriction)
            factors += [trans]
            trans_set.add(trans)

    monomial = format["multiply"](factors)
#    coeff_qe = format["multiply"](factors_qe)
#    if coeff_qe: coeff_qe = [coeff_qe]
#    else: coeff_qe = []

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
                    acs = format["add"]([str(c.index([], mapped_a, [], b)), str(offset)])
                    try: acs = "%d" % eval(acs)
                    except: pass
                    coefficient = format["coeff"] + format["matrix access"](str(c.n1.index), acs)
                else:
                    coefficient = format["coeff"] + format["matrix access"](c.n1.index, c.index([], mapped_a, [], b))
                for l in range(len(c.ops)):
                    op = c.ops[len(c.ops) - 1 - l]
                    if op == Operators.INVERSE:
                        coefficient = format["inverse"](coefficient)
                    elif op == Operators.MODULUS:
                        coefficient = format["absolute value"](coefficient)
                    elif op == Operators.SQRT:
                        coefficient = format["sqrt"](coefficient)
                factors += [coefficient]
        for t in tensor.transforms:
            if t.index0.type == Index.AUXILIARY_G or t.index1.type == Index.AUXILIARY_G:
#                print "t.type: ", t.type
#                print "t.restriction: ", t.restriction
#                print "t.index0: ", t.index0([],a,[],b)
#                print "t.index1: ", t.index1([],a,[],b)
                trans = format["transform"](t.type, t.index0([], a, [], b), \
                                                t.index1([], a, [], b), t.restriction)
                factors += [trans]
                trans_set.add(trans)

        terms += [format["multiply"](factors)]

    f_in = format["add"](terms)
    if f_in: f_in = [format["grouping"](f_in)]
    else: f_in = []

    return (f_out, f_in, trans_set)

def generate_terms(num_ips, tensor_numbers, tensors, format, psi_name_map, weight_name_map, opt_level):
    " Generate terms, returns a dictionary and the set of transformations used"

    list_indices    = format["free secondary indices"]
    primary_indices = [format["first free index"], format["second free index"]]
    format_multiply = format["multiply"]
    format_group    = format["grouping"]
    format_add      = format["add"]
    format_new_line = format["new line"]

    terms = {}

#    geo_num = len(geo_terms.keys()) + 1
    trans_set = Set()

    inv_w_name_map = {}
    # Invert weight name map
    for key in weight_name_map:
        vals = weight_name_map[key]
        if vals:
            for key_val in vals:
                inv_w_name_map[key_val] = key

    # Loop tensors and generate terms
    for tensor_number in tensor_numbers:
#            print "tensor number: ", tensor_number
        tensor = tensors[tensor_number]

        term = Term(tensor, list_indices)

        # Get rank and dims of primary indices
        irank, idims = tensor.i.rank, tensor.i.dims

        # Get monomial and compute macro dimensions in case of restricted basisfunctions
        monomial = tensor.monomial
        macro_idims = compute_macro_idims(monomial, idims, irank)

        # Get Psi indices, list of primary and secondary indices e.g. [[i0, a0], [i1, a1]]
        indices = term.indices
        vindices = term.vindices

        # Get secondary and auxiliary multiindices
        aindices, b0indices, bgindices = term.reduced_aindices, tensor.b0.indices, tensor.bg.indices

        # Get weight, take into account that we might only have 1 IP
        weight = format["weight"](tensor_number)
        if tensor_number in inv_w_name_map:
            weight = format["weight"](inv_w_name_map[tensor_number])
        if num_ips > 1:
            weight += format["array access"](format["integration points"])

        prim_dims = term.prim_dims
        prim_vars = term.prim_vars
        prim_dofs = term.prim_dofs
        sec_indices = term.secondary_indices

        for a in aindices:
            matrix_entry = {}
            # Change secondary index value to loop indices, for basis function indices
            for i in range(len(sec_indices)):
                a[sec_indices[i].index] = term.secondary_loops[i][0]
            for b0 in b0indices:

                a_map = {}
                sec_loop = [[],[]]
                R = []
                for psi_indices in indices:
#                        print "a: ", a
#                        print "psi indices: ", psi_indices
                    dof_map, dof_range, entry = generate_psi_entry2(num_ips, tensor_number, a,\
                      b0, psi_indices, vindices, psi_name_map, opt_level, format, tensor.qei)

#                    print "dof_map: ", dof_map
#                    print "dof_range: ", dof_range
#                    print "entry: ", entry

                    # If entry is '0' then the value is zero always
                    if entry:
                        R.append(entry)
                    if dof_map[0] in primary_indices:
#                            print "prim index: ", primary_indices.index(dof_map[0])
                        index_num = primary_indices.index(dof_map[0])
                        # All non-zero dofs
                        if dof_range == -1:
                            prim_dims[index_num] = idims[index_num]
                        else:
                            prim_dims[index_num] = dof_range
                        prim_vars[index_num] = dof_map[0]
                        prim_dofs[index_num] = dof_map[1]
                    else:
                        for li in term.secondary_loops:
                            # Do not add if dof_range == 1
                            if dof_map[0] == li[0]:
                                # All non-zero dofs
#                                    print "dof_map: ", dof_map
                                if dof_range == -1:
                                    sec_loop[0] += [li[2]]
                                    sec_loop[1] += [dof_map[0]]
                                elif not dof_range == 1:
                                    sec_loop[0] += [dof_range]
                                    sec_loop[1] += [dof_map[0]]

                    if not dof_map[0] in a_map:
                        a_map[dof_map[0]] = dof_map[1]
                    elif not a_map[dof_map[0]] == dof_map[1]:
                        raise RuntimeError, "Something is very wrong with secondary index map"

#                    print "dofs: ", prim_dofs

                # Generate entry name  for element tensor
                # FIXME: quadrature only support Functionals and Linear and Bilinear forms
                name = ""
                if (irank == 0):
                    # Entry is zero because functional is a scalar value
                    entry = "0"
                    # Generate name
                    name =  entry
                elif (irank == 1):
                    # Generate entry
                    entry = ""
                    for i in range(irank):
                        for v in monomial.basisfunctions:
                            if v.index.type == Index.PRIMARY and v.index.index == i:
                                if v.restriction == Restriction.MINUS:
                                    entry = format_add([prim_dofs[i], str(idims[0])])
                                else:
                                    entry = prim_dofs[i]
                                break

                    # Generate name
                    name = entry
                elif (irank == 2):
                    entry = []
                    # Generate entry
                    for i in range(irank):
                        for v in monomial.basisfunctions:
                            if v.index.type == Index.PRIMARY and v.index.index == i:
                                if v.restriction == Restriction.MINUS:
                                    entry += [format["grouping"](format_add([prim_dofs[i], str(idims[i])]))]
                                else:
                                    entry += [prim_dofs[i]]
                                break
#                        print "entry: ", entry
                    entry[0] = format_multiply([entry[0], str(macro_idims[1])])
                    name =  format_add(entry)
                else:
                    raise RuntimeError, "Quadrature only support Functionals and Linear and Bilinear forms"
#                    p_dof = tuple(prim_dofs)
#                print "name: ", name

                mapped_a = [index for index in a]
#                print "mapped_a: ", mapped_a
#                print "a_map: ", a_map
                # Map secondary indices to non-zero mapped indices
                for i in range(len(mapped_a)):
                    a_index = mapped_a[i]
                    if a_index in a_map:
                        mapped_a[i] = a_map[a_index]
#                print "a: ", mapped_a
                # Get geometry terms from inside sum, and outside sum
                geo_out, geo_in, t_set = generate_factor(tensor, a, mapped_a, bgindices, b0, format, opt_level)

                if 1 < len(geo_in):
                    geo_in = [format_group(format_add(geo_in))]

                d = ""
                if tensor.determinants:
                    d0 = [format["power"](format["determinant"](det.restriction),
                                      det.power) for det in tensor.determinants]
                    d = [format["multiply"]([format["scale factor"]] + d0)]
                else:
                    d = [format["scale factor"]]

                geo = format_multiply(geo_out + geo_in + d)

                val = ""
                if R:
                    val = format_multiply(R + [weight])
                else:
                    val = weight

                # We only update the sets if we're going to add this block
                # of code i.e. None is not in R
                trans_set.add(d[0])
                trans_set = trans_set | t_set

                # Add contribution to code
                prim_loop = (tuple(prim_dims), tuple(prim_vars))
                sec_loop = tuple([tuple(s) for s in sec_loop])
                if (prim_loop, sec_loop) not in terms:
                    terms[(prim_loop, sec_loop)] = {}
                    terms[(prim_loop, sec_loop)][name] = [format_multiply([val, geo])]
                else:
                    if name not in terms[(prim_loop, sec_loop)]:
                        terms[(prim_loop, sec_loop)][name] = [format_multiply([val, geo])]
                    else:
                        terms[(prim_loop, sec_loop)][name].append( format_multiply([val, geo]) )

#            print non_zero_columns, cols_name_map
#            print inv_w_name_map

    return (terms, trans_set)










