"Code generator for quadrature representation"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2007-03-16 -- 2007-05-07"
__copyright__ = "Copyright (C) 2004-2007 Kristian B. Oelgaard"
__license__  = "GNU GPL Version 2"

# Modified by Anders Logg 2007

# Python modules
import numpy

# FFC common modules
from ffc.common.constants import *
from ffc.common.utils import *

# FFC language modules
from ffc.compiler.language.index import *
from ffc.compiler.language.restriction import *

# FFC code generation modules
from ffc.compiler.codegeneration.common.codegenerator import *
from ffc.compiler.codegeneration.common.utils import *

# FFC tensor representation modules
from ffc.compiler.representation.tensor.multiindex import *

# Should be in dictionary!!
index_names = {0: lambda i: "f%s" %(i), 1: lambda i: "p%s" %(i), 2: lambda i: "s%s" %(i),\
               4: lambda i: "fu%s" %(i), 5: lambda i: "pj%s" %(i), 6: lambda i: "c%s" %(i), 7: lambda i: "a%s" %(i)}

# Should be in dictionary!!
#index_names = {0: lambda i: "fix_%s" %(i), 1: lambda i: "prim_%s" %(i), 2: lambda i: "sec_%s" %(i),\
#               3: lambda i: "aux_%s" %(i), 4: lambda i: "func_%s" %(i), 5: lambda i: "proj_%s" %(i),\
#               6: lambda i: "cons_%s" %(i), 7: lambda i: "aux0_%s" %(i)}


class QuadratureGenerator(CodeGenerator):
    "Code generator for for tensor representation"

    def __init__(self):
        "Constructor"

        # Initialize common code generator
        CodeGenerator.__init__(self)

    def generate_cell_integral(self, form_representation, sub_domain, format):
        """Generate dictionary of code for cell integral from the given
        form representation according to the given format"""

        # Object to control the code indentation
        Indent = IndentControl()

        # Extract tensors
        tensors = form_representation.cell_tensor
        if len(tensors) == 0:
            return None

        # Generate code for element tensor(s)
        code = []
        code += [Indent.indent(format["comment"]("Compute element tensor"))]
        code += self.__generate_element_tensor(tensors, Indent, format)

        return {"tabulate_tensor": code}

    def generate_exterior_facet_integral(self, form_representation, sub_domain, format):
        """Generate dictionary of code for exterior facet integral from the given
        form representation according to the given format"""

        # Object to control the code indentation
        Indent = IndentControl()

        # Prefetch formats to speed up code generation
        format_block_begin  = format["block begin"]
        format_block_end    = format["block end"]

        # Extract tensors
        tensors = form_representation.exterior_facet_tensors
        if len(tensors) == 0:
            return None

        # Generate code for element tensor(s)
        common = [""] + [format["comment"]("Compute element tensor for all facets")]

        num_facets = len(tensors)
        cases = [None for i in range(num_facets)]
        for i in range(num_facets):
            case = [format_block_begin]
            case += self.__generate_element_tensor(tensors[i], Indent, format)
            case += [format_block_end]
            cases[i] = case

        return {"tabulate_tensor": (common, cases)}
    
    def generate_interior_facet_integral(self, form_representation, sub_domain, format):
        """Generate dictionary of code for interior facet integral from the given
        form representation according to the given format"""

        # Object to control the code indentation
        Indent = IndentControl()

        # Prefetch formats to speed up code generation
        format_block_begin  = format["block begin"]
        format_block_end    = format["block end"]

        # Extract tensors
        tensors = form_representation.interior_facet_tensors
        if len(tensors) == 0:
            return None

        # Generate code for element tensor(s)
        common = [""] + [format["comment"]("Compute element tensor for all facet-facet combinations")]
        num_facets = len(tensors)
        cases = [[None for j in range(num_facets)] for i in range(num_facets)]
        for i in range(num_facets):
            for j in range(num_facets):
                case = [format_block_begin]
                case += self.__generate_element_tensor(tensors[i][j], Indent, format)
                case += [format_block_end]
                cases[i][j] = case

        return {"tabulate_tensor": (common, cases)}

    def __generate_element_tensor(self, tensors, Indent, format):
        "Construct quadrature code for element tensors"

        # Initialize code segments
        tabulate_code = []
        element_code = []

        # Prefetch formats to speed up code generation
        format_comment      = format["comment"]
        format_loop         = format["loop"]
        format_ip           = format["integration points"]
        format_block_begin  = format["block begin"]
        format_block_end    = format["block end"]

        group_tensors = self.__equal_num_quadrature_points(tensors)
#        print "group_tensors: ", group_tensors
        # Number of tensors to evaluate
        num_tensors = len(tensors)
#        print "num_tensors: ", num_tensors

        for points in group_tensors:
            for i in group_tensors[points]:
        # Loop tensors to generate tables
#        for i in range(num_tensors):

            # Tabulate variables:
            # Tabulate the quadrature weights
                tabulate_code += self.__tabulate_weights(tensors[i].quadrature.weights, i, Indent, format)

            # Tabulate values of basis functions and their derivatives at quadrature points
                tabulate_code += self.__tabulate_psis(tensors[i], i, Indent, format)

        # Reset values of the element tensor (assuming same dimensions for all tensors)
            tabulate_code += self.__reset_element_tensor(tensors[0], Indent, format)

        # Loop tensors to generate quadrature loops
#        for i in range(num_tensors):
            # Loop all quadrature points
            element_code += [Indent.indent(format_comment("Loop quadrature points (tensor/monomial term %d)" %(i,)))]
            element_code += [Indent.indent(format_loop(format_ip, 0, len(tensors[i].quadrature.weights)))]
            element_code += [Indent.indent(format_block_begin)]

            # Increase indentation
            Indent.increase()

            # Generate the element tensor
            for i in group_tensors[points]:
                element_code += self.__element_tensor(tensors[i], i, Indent, format)

            # Decrease indentation
            Indent.decrease()
            # End the quadrature loop
            element_code += [Indent.indent(format_block_end)]

            if i + 1 < len(tensors):
                element_code += [""]

        return tabulate_code + element_code

    def __tabulate_weights(self, weights, tensor_number, Indent, format):
        "Generate table of quadrature weights"

        code = []
        
        # Prefetch formats to speed up code generation
        format_floating_point = format["floating point"]

        # Get number of weights
        num_weights = len(weights)

        code += [Indent.indent(format["comment"]\
                ("Array of quadrature weights (tensor/monomial term %d)" %(tensor_number,) ))]

        # Create variable name
        name = format["table declaration"] + format["weights"](tensor_number, str(num_weights))
        value = format["block"](format["separator"].join([format_floating_point(weights[i])\
                 for i in range(num_weights)]))

        code += [(Indent.indent(name), value)]

        return code + [""]

    def __tabulate_psis(self, tensor, tensor_number, Indent, format):
        "Tabulate values of basis functions and their derivatives at quadrature points"

        code = []

        format_psis = format["psis"]
        format_floating_point = format["floating point"]
        format_block_begin = format["block begin"]
        format_block_end = format["block end"]

        # Get list of psis
        psis = tensor.Psis

        code += [Indent.indent(format["comment"]\
                ("Values of shapefunctions and their derivatives at quadrature points"))]
        code += [Indent.indent(format["comment"]\
                ("Format: [quadrature points][dofs] (tensor/monomial term %d)" % (tensor_number,)))]

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

                    (name, multi_index) = self.__generate_psi_declaration(tensor_number, indices, vindex,\
                                a, b, num_quadrature_points, num_dofs, format)

                    names += [name]
                    multi_indices += [multi_index]

            # Remove redundant names and entries in the psi table
            names, multi_indices = self.__extract_unique(names, multi_indices)

            # Loop names and tabulate psis
            for i in range(len(names)):
                # Get values from psi tensor, should have format values[dofs][quad_points]
                vals = values[tuple(multi_indices[i])]

                # Generate array of values (FIAT returns [dof, quad_points] transpose to [quad_points, dof])
                value = tabulate_matrix(numpy.transpose(vals), format)

                code += [(Indent.indent(names[i]), Indent.indent(value))]# + [""]

        return code

    def __reset_element_tensor(self, tensor, Indent, format):
        "Reset the entries of the local element tensor"

        code = []

        # Get rank and dims of primary indices
        irank = tensor.i.rank
        idims = tensor.i.dims

        # Create macro dimensions in case of restricted basisfunctions
        macro_idims = [dim for dim in idims]
        monomial = tensor.monomial
        for i in range(irank):
            for v in monomial.basisfunctions:
                if v.index.type == Index.PRIMARY and v.index.index == i:
                    if v.restriction != None:
                        macro_idims[i] = idims[i] * 2
                        break

        # Generate value
        value = format["floating point"](0.0)

        # FIXME: quadrature only supports Linear and Bilinear forms
        if (irank == 1):
            code += [Indent.indent(format["comment"]\
                    ("Reset values of the element tensor block"))]

            # Generate entry and name
            entry = format["first free index"]
            name =  format["element tensor quad"] + format["array access"](entry)

            # Create boundaries for loop
            boundaries = [0, macro_idims[0]]
            code += self.__generate_loop(name, value, boundaries, Indent, format)

        elif (irank == 2):
            code += [Indent.indent(format["comment"]\
                    ("Reset values of the element tensor block"))]

            # Generate entry and name
            entry = format["add"]([format["multiply"]([format["first free index"], "%d" %macro_idims[1]]),\
                                                       format["second free index"]])
            name =  format["element tensor quad"] + format["array access"](entry)

            # Create boundaries for loop
            boundaries = [0, macro_idims[0], 0, macro_idims[1]]
            code += self.__generate_loop(name, value, boundaries, Indent, format)
        else:
            raise RuntimeError, "Quadrature only supports Linear and Bilinear forms"

        return code + [""]

    def __element_tensor(self, tensor, tensor_number, Indent, format):
        "Generate loop over primary indices"

        code = []

        # Prefetch formats to speed up code generation
        format_add      = format["add"]
        format_multiply = format["multiply"]
        format_group    = format["grouping"]
        format_new_line = format["new line"]

        dic_indices = {0:format["first free index"], 1:format["second free index"]}

        # Get rank and dims of primary indices
        irank, idims = tensor.i.rank, tensor.i.dims

        # Initialise macro dimensions and get monomial
        macro_idims = [dim for dim in idims]
        monomial = tensor.monomial
        for i in range(irank):
            for v in monomial.basisfunctions:
                if v.index.type == Index.PRIMARY and v.index.index == i:
                    if v.restriction != None:
                        macro_idims[i] = idims[i] * 2
                        break

        # Get Psi indices, list of primary and secondary indices e.g. [[i0, a0], [i1, a1]]
        indices = [psi[1] for psi in tensor.Psis]
        vindices = [psi[2] for psi in tensor.Psis]

        # Get secondary and auxiliary multiindices
        aindices, b0indices, bgindices = tensor.a.indices, tensor.b0.indices, tensor.bg.indices

        # Compute scaling
        weight = [format["weights"](tensor_number, format["integration points"])]

        # Generate brackets of geometry and reference terms (optimised)
        values = []
        for a in aindices:
            r, g = [], []
            for b0 in b0indices:
                r += [format_multiply([self.__generate_psi_entry(tensor_number, a,\
                                b0, psi_indices, vindices, format)\
                       for psi_indices in indices] + weight)]

            if 1 < len(r):
                ref = format_group(format_add(r))
            else:
                ref = r[0]

            for bg in bgindices:
                g += [format_multiply(self.__generate_factor(tensor, a, bg, format))]
            if 1 < len(g):
                geo = format_group(format_add(g))
            else:
                geo = g[0]
            values += [format_multiply([ref,geo]) + format_new_line]
        value = format_add(values)

        # Generate value (expand multiplication - not optimised)
#        values = []
#        for a in aindices:
#            for b0 in b0indices:
#                for bg in bgindices:
#                    factor = self.__generate_factor(tensor, a, bg, format)
#                    values += [format_multiply([self.__generate_psi_entry(tensor_number, a,\
#                                                          b0, psi_indices, vindices, format)\
#                         for psi_indices in indices] + weight + factor) + format["new line"]]

#        value = format_add(values)

        # FIXME: quadrature only supports Linear and Bilinear forms
        if (irank == 1):

            # Generate entry
            for i in range(irank):
                for v in monomial.basisfunctions:
                    if v.index.type == Index.PRIMARY and v.index.index == i:
                        if v.restriction == Restriction.MINUS:
                            entry = format_add([dic_indices[i], str(idims[0])])
                        else:
                            entry = dic_indices[i]
                        break

            # Generate name
            name =  format["element tensor quad"] + format["array access"](entry)
            code += [Indent.indent(format["comment"]("Compute block entries (tensor/monomial term %d)" % (tensor_number,)))]

            # Create boundaries for loop
            boundaries = [0, idims[0]]
            code += self.__generate_loop(name, value, boundaries, Indent, format, format["add equal"])

        elif (irank == 2):

            entry = []
            # Generate entry
            for i in range(irank):
                for v in monomial.basisfunctions:
                    if v.index.type == Index.PRIMARY and v.index.index == i:
                        if v.restriction == Restriction.MINUS:
                            entry += [format["grouping"](format_add([dic_indices[i], str(idims[i])]))]
                        else:
                            entry += [dic_indices[i]]
                        break

            entry[0] = format_multiply([entry[0], str(macro_idims[1])])
            name =  format["element tensor quad"] + format["array access"](format_add(entry))

            code += [Indent.indent(format["comment"]("Compute block entries (tensor/monomial term %d)" % (tensor_number,)))]

            # Create boundaries for loop
            boundaries = [0, idims[0], 0, idims[1]]
            code += self.__generate_loop(name, value, boundaries, Indent, format, format["add equal"])
        else:
            raise RuntimeError, "Quadrature only supports Linear and Bilinear forms"

        return code

    def __generate_psi_declaration(self, tensor_number, psi_indices, vindex, aindices, bindices,\
                                         num_quadrature_points, num_dofs, format):

        # Prefetch formats to speed up code generation
        format_secondary_index  = format["secondary index"]

        # Should be in dictionary!!
#        index_names = {0: lambda i: "fix_%s" %(i), 1: lambda i: "prim_%s" %(i), 2: lambda i: "sec_%s" %(i),\
#                 3: lambda i: "aux_%s" %(i), 4: lambda i: "func_%s" %(i), 5: lambda i: "proj_%s" %(i),\
#                 6: lambda i: "cons_%s" %(i), 7: lambda i: "aux0_%s" %(i)}

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

    def __generate_psi_entry(self, tensor_number, aindices, bindices, psi_indices, vindices, format):

        # Prefetch formats to speed up code generation
        format_secondary_index  = format["secondary index"]

        # Should be in dictionary!!
#        index_names = {0: lambda i: "fix_%s" %(i), 1: lambda i: "prim_%s" %(i), 2: lambda i: "sec_%s" %(i),\
#                 3: lambda i: "aux_%s" %(i), 4: lambda i: "func_%s" %(i), 5: lambda i: "proj_%s" %(i),\
#                 6: lambda i: "cons_%s" %(i), 7: lambda i: "aux0_%s" %(i)}


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

    def __generate_factor(self, tensor, a, b, format):
        "Generate code for the value of entry a of geometry tensor G"

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
#        monomial = format["multiply"](factors)
#        if monomial: f0 = [monomial]
#        else: f0 = []
    
        # Compute sum of monomials inside sum
#        terms = []
#        for b in G.b.indices:
#            factors = []
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
# End from tensorgenerator

        if tensor.determinant:
            d0 = format["power"](format["determinant"], tensor.determinant)
            d = format["multiply"]([format["scale factor"], d0])
        else:
            d = format["scale factor"]

        factors += [d]

        return factors

    def __generate_loop(self, name, value, boundaries, Indent, format, connect = None):
        "This function generates a loop over a vector or matrix."

        code = []

        # Prefetch formats to speed up code generation
        format_loop     = format["loop"]

        if (len(boundaries) == 2):
            # Loop index
            code += [Indent.indent(format_loop(format["first free index"], boundaries[0], boundaries[1]))]
            # Increase indentation
            Indent.increase()

            if connect:
                code += [connect(Indent.indent(name), value)]
            else:
                code += [(Indent.indent(name), value)]

            # Decrease indentation
            Indent.decrease()

        elif (len(boundaries) == 4):
            # Loop first primary index
            code += [Indent.indent(format_loop(format["first free index"], boundaries[0], boundaries[1]))]
            code += [Indent.indent(format["block begin"])]

            # Increase indentation
            Indent.increase()

            # Loop second primary index
            code += [Indent.indent(format_loop(format["second free index"], boundaries[2], boundaries[3]))]

            # Increase indentation
            Indent.increase()

            if connect:
                code += [connect(Indent.indent(name), value)]
            else:
                code += [(Indent.indent(name), value)]

            # Decrease indentation
            Indent.decrease()
            # Decrease indentation
            Indent.decrease()
            code += [Indent.indent(format["block end"])]

        else:
            raise RuntimeError, "This function can only generate a loop for a vector or a matrix"

        return code

    def __extract_unique(self, aa, bb):
        "Remove redundant names and entries in the psi table"

        uaa = []
        ubb = []

        for i in range(len(aa)):
            a = aa[i]
            if not a in uaa:
                uaa += [a]
                ubb += [bb[i]]
        return (uaa, ubb)

    def __equal_num_quadrature_points(self, tensors):

        group_tensors = {}
        for i in range(len(tensors)):
            tens = tensors[i]
            num_points = len(tens.quadrature.weights)
            if num_points in group_tensors:
                group_tensors[num_points] += [i]
            else:
                group_tensors[num_points] = [i]

        return group_tensors
