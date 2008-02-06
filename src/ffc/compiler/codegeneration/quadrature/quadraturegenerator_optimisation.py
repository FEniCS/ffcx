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

def generate_factor_old(tensor, a, b, format):
    "Optimise level 0 and 1"
    # Compute product of factors outside sum
    factors = []
    trans_set = Set()

    for j in range(len(tensor.coefficients)):
        c = tensor.coefficients[j]
        if not c.index.type == Index.AUXILIARY_G:
            offset = tensor.coefficient_offsets[c]
            coefficient = format["coefficient"](c.n1.index, c.index([], a, [], [])+offset)
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
                                                    t.index1([], a, [], []), \
                                                    t.restriction)
            factors += [trans]
            trans_set.add(trans)

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
                elif op == Operators.MODULUS:
                    coefficient = format["absolute value"](coefficient)
                elif op == Operators.SQRT:
                    coefficient = format["sqrt"](coefficient)
            factors += [coefficient]
    for t in tensor.transforms:
        if t.index0.type == Index.AUXILIARY_G or t.index1.type == Index.AUXILIARY_G:
            trans= format["transform"](t.type, t.index0([], a, [], b), \
                                                        t.index1([], a, [], b), \
                                                        t.restriction)
            factors += [trans]
            trans_set.add(trans)

    if tensor.determinants:
        d0 = [format["power"](format["determinant"](det.restriction),
                                          det.power) for det in tensor.determinants]
        d = format["multiply"]([format["scale factor"]] + d0)
    else:
        d = format["scale factor"]

    trans_set.add(d)

    return ([format["multiply"](factors + [d])], trans_set)

def generate_factor(tensor, a, bgindices, format):
    "Optimise level 2"

    trans_set = Set()

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
                elif op == Operators.MODULUS:
                    coefficient = format["MODULUSolute value"](coefficient)
                elif op == Operators.SQRT:
                    coefficient = format["sqrt"](coefficient)
            factors += [coefficient]
    for t in tensor.transforms:
        if not (t.index0.type == Index.AUXILIARY_G or  t.index1.type == Index.AUXILIARY_G):
            trans = format["transform"](t.type, t.index0([], a, [], []), \
                                                    t.index1([], a, [], []), \
                                                    t.restriction)
            factors += [trans]
            trans_set.add(trans)

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
                    elif op == Operators.MODULUS:
                        coefficient = format["absolute value"](coefficient)
                    elif op == Operators.SQRT:
                        coefficient = format["sqrt"](coefficient)
                factors += [coefficient]
        for t in tensor.transforms:
            if t.index0.type == Index.AUXILIARY_G or t.index1.type == Index.AUXILIARY_G:
                trans = format["transform"](t.type, t.index0([], a, [], b), \
                                                        t.index1([], a, [], b), \
                                                        t.restriction)
                factors += [trans]
                trans_set.add(trans)

        terms += [format["multiply"](factors)]

    f_in = format["add"](terms)
    if f_in: f_in = [format["grouping"](f_in)]
    else: f_in = []

    return (f_out, f_in, trans_set)

def generate_factor3(tensor, a, bgindices, format):
    "Optimise level 3 and 4"

    # Compute product of factors outside sum
    trans_set = Set()
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
                    elif op == Operators.MODULUS:
                        coefficient = format["absolute value"](coefficient)
                    elif op == Operators.SQRT:
                        coefficient = format["sqrt"](coefficient)
                factors += [coefficient]
        for t in tensor.transforms:
            if t.index0.type == Index.AUXILIARY_G or t.index1.type == Index.AUXILIARY_G:
                trans = format["transform"](t.type, t.index0([], a, [], b), \
                                                t.index1([], a, [], b), t.restriction)
                factors += [trans]
                trans_set.add(trans)

        terms += [format["multiply"](factors)]

    f_in = format["add"](terms)
    if f_in: f_in = [format["grouping"](f_in)]
    else: f_in = []

    return (f_out, f_in, trans_set)

def generate_factor5(tensor, a, bgindices, b0, format):
    "Optimise level 5"

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
            access = str(c.index([], a, [], []))
            if c.index in tensor.qei:
                for v in tensor.monomial.basisfunctions:
                    if v.index == c.index:
#                        print "v.component: ", v.component
                        if v.component:
                            if len(v.component) == 1:
                                access = qe_access[v.component[0]([], a, b0, [])]
                            else:
                                raise RuntimeError("Error, more than one component index!!")
                        else:
                            access = qe_access[0]
#                            raise RuntimeError("Error!!")

            if offset:
                coefficient = format["coeff"] + format["matrix access"](str(c.n1.index),\
                              format["add"]([access, str(offset)]))
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
    coeff_qe = format["multiply"](factors_qe)
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
                    coefficient = format["coeff"] + format["matrix access"](str(c.n1.index),\
                                  format["add"]([str(c.index([], a, [], b)), str(offset)]))
                else:
                    coefficient = format["coeff"] + format["matrix access"](c.n1.index, c.index([], a, [], b))
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
                trans = format["transform"](t.type, t.index0([], a, [], b), \
                                                t.index1([], a, [], b), t.restriction)
                factors += [trans]
                trans_set.add(trans)

        terms += [format["multiply"](factors)]

    f_in = format["add"](terms)
    if f_in: f_in = [format["grouping"](f_in)]
    else: f_in = []

    return (f_out, f_in, trans_set)

def values_level_0(indices, vindices, aindices, b0indices, bgindices, tensor, tensor_number, weight, format, name_map):
    # Generate value (expand multiplication - optimise level 0)
    format_multiply = format["multiply"]
    format_new_line = format["new line"]

    values = []
    trans_set = Set()

    for a in aindices:
        for b0 in b0indices:
            for bg in bgindices:
                factor, t_set = generate_factor_old(tensor, a, bg, format)
                trans_set = trans_set | t_set

                values += [format_multiply([generate_psi_entry(tensor_number, a,\
                           b0, psi_indices, vindices, name_map, format) for psi_indices in indices] +\
                           weight + factor) + format_new_line]

    return (values, [], trans_set)

def values_level_1(indices, vindices, aindices, b0indices, bgindices, tensor, tensor_number, weight, format, name_map):
    # Generate brackets of geometry and reference terms (terms outside sum are
    # multiplied with each term in the sum  - optimise level 1)
    format_multiply = format["multiply"]
    format_group = format["grouping"]
    format_add = format["add"]
    format_new_line = format["new line"]

    values = []
    trans_set = Set()

    for a in aindices:
        r = []
        for b0 in b0indices:
            r += [format_multiply([generate_psi_entry(tensor_number, a,\
                  b0, psi_indices, vindices, name_map, format) for psi_indices in indices] + weight)]

        if 1 < len(r):
            ref = format_group(format_add(r))
        else:
            ref = r[0]

        geo = []
        for bg in bgindices:
            factor, t_set = generate_factor_old(tensor, a, bg, format)
            trans_set = trans_set | t_set
            geo += [format_multiply(factor)]

        if 1 < len(geo):
            geo = format_group(format_add(geo))
        else:
            geo = geo[0]
        values += [format_multiply([ref,geo]) + format_new_line]

    return (values, [], trans_set)

def values_level_2(indices, vindices, aindices, b0indices, bgindices, tensor, tensor_number, weight, format, name_map):
    # Generate brackets of geometry and reference terms, distinguish between terms outside and 
    # inside sum (geometry only). - optimise level 2)

    format_multiply = format["multiply"]
    format_group = format["grouping"]
    format_add = format["add"]
    format_new_line = format["new line"]

    values = []
    trans_set = Set()

    for a in aindices:
        r = []
        for b0 in b0indices:
            r += [format_multiply([generate_psi_entry(tensor_number, a,\
                  b0, psi_indices, vindices, name_map, format) for psi_indices in indices] + weight)]

        if 1 < len(r):
            ref = format_group(format_add(r))
        else:
            ref = r[0]

        # Get geometry terms from inside sum, and outside sum
        geo_out, geo_in, t_set = generate_factor(tensor, a, bgindices, format)
        trans_set = trans_set | t_set

        if 1 < len(geo_in):
            geo_in = [format_group(format_add(geo_in))]

        if tensor.determinants:
            d0 = [format["power"](format["determinant"](det.restriction),
                                  det.power) for det in tensor.determinants]
            d = [format["multiply"]([format["scale factor"]] + d0)]
        else:
            d = [format["scale factor"]]

        trans_set.add(d[0])

        geo = format_multiply(geo_out + geo_in + d)
        values += [format_multiply([ref,geo]) + format_new_line]

    return (values, [], trans_set)

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
    trans_set = Set()
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
                  b0, psi_indices, vindices, name_map, format) for psi_indices in indices] + weight)]

        if 1 < len(r):
            ref = format_group(format_add(r))
        else:
            ref = r[0]

        # Get geometry terms from inside sum, and outside sum
        geo_out, geo_in, t_set = generate_factor3(tensor, a, bgindices, format)
        trans_set = trans_set | t_set
        if 1 < len(geo_in):
            geo_in = [format_group(format_add(geo_in))]

        if tensor.determinants:
            d0 = [format["power"](format["determinant"](det.restriction),
                                  det.power) for det in tensor.determinants]
            d = [format["multiply"]([format["scale factor"]] + d0)]
        else:
            d = [format["scale factor"]]

        trans_set.add(d[0])

        geo = format_multiply(geo_out + geo_in + d)

        vals += [format_multiply([ref,geo]) + format_new_line]

    # Only use values that are unique
    for val in vals:
        if not val in values:
            values += [val]
#    print "\nvalues[0]: ", values[0]

    return (values, secondary_loop, trans_set)

def values_level_4(indices, vindices, aindices, b0indices, bgindices, tensor, tensor_number, weight, format, name_map):
    # Generate value - optimise level 3 (Based on level 2 but use loops to multiply reference and geo tensor)

    format_multiply = format["multiply"]
    format_group    = format["grouping"]
    format_add      = format["add"]
    format_new_line = format["new line"]
    format_space    = format["space"]
    format_psis     = format["psis"]
    format_if       = format["if"]
    format_and      = format["logical and"]
    format_abs      = format["absolute value"]
    format_greater  = format["greater than"]
    format_eps      = format["epsilon"]
    format_float    = format["floating point"]

    # Get list of free secondary loop indices
    list_indices = format["free secondary indices"]
    values = []
    vals = []
    secondary_loop = []
    sec_indices = []
    trans_set = Set()
    num_if = 0
    if_statements = []
    if_s = []
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

            # Pick index and generate information for loop generation
            index_name = old_ind[index.index]
            secondary_loop += [[index_name, 0, len(index.range)]]

    for a in aindices:
        # Change secondary index value to loop indices, for basis function indices
        for i in range(len(sec_indices)):
            a[sec_indices[i].index] = secondary_loop[i][0]
        ref = []
        for b0 in b0indices:
            entries = [generate_psi_entry(tensor_number, a,\
                  b0, psi_indices, vindices, name_map, format) for psi_indices in indices]
            if_s.append( format_if + format_group(format_and.join([format_abs(e) +\
                                    format_greater + format_float(format_eps) for e in entries])) )
            ref += [format_multiply(entries + weight)]

#            ref += [format_multiply([generate_psi_entry(tensor_number, a,\
#                  b0, psi_indices, vindices, name_map, format) for psi_indices in indices] + weight)]

        # Get geometry terms from inside sum, and outside sum
        geo_out, geo_in, t_set = generate_factor3(tensor, a, bgindices, format)
        trans_set = trans_set | t_set
        if 1 < len(geo_in):
            geo_in = [format_group(format_add(geo_in))]

        if tensor.determinants:
            d0 = [format["power"](format["determinant"](det.restriction),
                                  det.power) for det in tensor.determinants]
            d = [format["multiply"]([format["scale factor"]] + d0)]
        else:
            d = [format["scale factor"]]

        trans_set.add(d[0])

        geo = format_multiply(geo_out + geo_in + d)

        vals += [format_multiply([re,geo]) + format_new_line for re in ref]

    # Only use values that are unique
    for i in range(len(vals)):
        val = vals[i]
        if_ = if_s[i]
        if not val in values:
            values.append(val)
            if_statements.append(if_)

    return (values, secondary_loop, trans_set, if_statements)

def values_level_5(indices, vindices, aindices, b0indices, bgindices, tensor, tensor_number, weight, format, name_map):
    # Generate value - optimise level 3 (Based on level 2 but use loops to multiply reference and geo tensor)

    format_multiply = format["multiply"]
    format_group    = format["grouping"]
    format_add      = format["add"]
    format_new_line = format["new line"]
    format_space    = format["space"]
    format_psis     = format["psis"]
    format_if       = format["if"]
    format_and      = format["logical and"]
    format_abs      = format["absolute value"]
    format_greater  = format["greater than"]
    format_eps      = format["epsilon"]
    format_float    = format["floating point"]

    # Get list of free secondary loop indices
    list_indices = format["free secondary indices"]

    values = []
    vals = []
    secondary_loop = []
    sec_indices = []
    trans_set = Set()
    if_statements = []
    if_s = []

    for index in vindices:
#        print "index: ", index
        if index.type == Index.SECONDARY and len(index.range) > 1 and not index in tensor.qei:
#        if index.type == Index.SECONDARY and len(index.range) > 1:
            sec_indices += [index]
            # Generate loop variables
            old_ind = [d for d in list_indices]
#            print "old ind: ", old_ind
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

    for a in aindices:
        # Change secondary index value to loop indices, for basis function indices
        for i in range(len(sec_indices)):
            a[sec_indices[i].index] = secondary_loop[i][0]
#        print "a: ", a
        ref = []
        for b0 in b0indices:
            entries = [generate_psi_entry(tensor_number, a,\
                  b0, psi_indices, vindices, name_map, format, tensor.qei) for psi_indices in indices]
#            entries = [generate_psi_entry(tensor_number, a,\
#                  b0, psi_indices, vindices, name_map, format) for psi_indices in indices]

            if_s.append( format_if + format_group(format_and.join([format_abs(e) +\
                                    format_greater + format_float(format_eps) for e in entries if e])) )
            if "" in entries:
                entries.remove("")
#            ref += [format_multiply(entries + weight)]
            re = format_multiply(entries + weight)

            # Get geometry terms from inside sum, and outside sum
            geo_out, geo_in, t_set = generate_factor5(tensor, a, bgindices, b0, format)
            trans_set = trans_set | t_set
            if 1 < len(geo_in):
                geo_in = [format_group(format_add(geo_in))]

            if tensor.determinants:
                d0 = [format["power"](format["determinant"](det.restriction),
                                  det.power) for det in tensor.determinants]
                d = [format["multiply"]([format["scale factor"]] + d0)]
            else:
                d = [format["scale factor"]]

            trans_set.add(d[0])
            geo = format_multiply(geo_out + geo_in + d)

#            vals += [format_multiply([re, geo]) + format_new_line for re in ref]
            vals.append(format_multiply([re, geo]) + format_new_line)

    # Only use values that are unique
    for i in range(len(vals)):
        val = vals[i]
        if_ = if_s[i]
        if not val in values:
            values.append(val)
            if_statements.append(if_)

    return (values, secondary_loop, trans_set, if_statements)

