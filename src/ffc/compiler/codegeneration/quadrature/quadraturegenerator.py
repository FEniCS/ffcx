"Code generator for quadrature representation"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2007-03-16 -- 2007-03-23"
__copyright__ = "Copyright (C) 2004-2007 Kristian B. Oelgaard"
__license__  = "GNU GPL Version 2"

# Python modules
import numpy

# FFC common modules
from ffc.common.constants import *
from ffc.common.utils import *

# FFC language modules
from ffc.compiler.language.index import *

# FFC code generation modules
from ffc.compiler.codegeneration.common.codegenerator import *

from ffc.compiler.representation.quadrature.multiindex import *

class QuadratureGenerator(CodeGenerator):
    "Code generator for for tensor representation"

    def __init__(self):
        "Constructor"

        # Initialize common code generator
        CodeGenerator.__init__(self)

        # Initialize the indentation size and increment
        self.indent_size = 0
        self.indent_increment = 2

    def generate_cell_integral(self, form_representation, sub_domain, format):
        """Generate dictionary of code for cell integral from the given
        form representation according to the given format"""

        # Extract tensors
        tensors = form_representation.cell_tensor
        if len(tensors) == 0:
            return None

        # Generate code for geometry tensor
#        code = self.__generate_geometry_tensors(tensors, format)
        
#        # Generate code for element tensor(s)
#        code += [""] + [format["comment"]("Compute element tensor")]
#        code += self.__generate_element_tensor(tensors, format)

        # Generate code for element tensor(s)
        code = []
        code += self.__generate_element_tensor(tensors, format)
#        for i in range(len(tensors)):
#            code += self.__generate_element_tensor(tensors[i], i, format)

        return {"tabulate_tensor": code}

    def generate_exterior_facet_integral(self, form_representation, sub_domain, format):
        """Generate dictionary of code for exterior facet integral from the given
        form representation according to the given format"""

        # Extract tensors
        tensors = form_representation.exterior_facet_tensors
#        if len(tensors) == 0:
#            return None

        # Generate code for geometry tensor (should be the same so pick first)
#        common = self.__generate_geometry_tensors(tensors[0], format)

        # Generate code for element tensor(s)
#        common += [""] + [format["comment"]("Compute element tensor for all facets")]
#        num_facets = len(tensors)
#        cases = [None for i in range(num_facets)]
#        for i in range(num_facets):
#            cases[i] = self.__generate_element_tensor(tensors[i], format)
        common = [""]
        cases = [""]
        return {"tabulate_tensor": (common, cases)}
    
    def generate_interior_facet_integral(self, form_representation, sub_domain, format):
        """Generate dictionary of code for interior facet integral from the given
        form representation according to the given format"""

        # Extract tensors
        tensors = form_representation.interior_facet_tensors
#        if len(tensors) == 0:
#            return None
        
        # Generate code for geometry tensor (should be the same so pick first)
#        common = self.__generate_geometry_tensors(tensors[0][0], format)

        # Generate code for element tensor(s)
#        common += [""] + [format["comment"]("Compute element tensor for all facet-facet combinations")]
#        num_facets = len(tensors)
#        cases = [[None for j in range(num_facets)] for i in range(num_facets)]
#        for i in range(num_facets):
#            for j in range(num_facets):
#                cases[i][j] = self.__generate_element_tensor(tensors[i][j], format)
        common = [""]
        cases = [""]
        return {"tabulate_tensor": (common, cases)}

    def __generate_element_tensor(self, tensors, format):
        "Construct quadrature code for element tensors"

        # Initialize code segments
        tabulate_code = []
        element_code = []

        format_comment = format["comment"]

        # Number of tensors to evaluate
        num_tensors = len(tensors)

        # Loop tensors to generate tables
        for i in range(num_tensors):

            # Tabulate variables,
            # Tabulate the derivatives of shapefunctions, used to generate the Jacobian
            # FIXME: Switched off until a nodemap is available (Jacobian)
#            tabulate_code += self.__tabulate_derivatives(tensors[i].Derivatives, i, format)

            # Tabulate the quadrature weights
            tabulate_code += self.__tabulate_weights(tensors[i].quadrature.weights, i, format)

            # Tabulate values of basis functions and their derivatives at quadrature points
            tabulate_code += self.__tabulate_psis(tensors[i], i, format)

        # Reset values of the element tensor (assuming same dimensions for the tensors)
        tabulate_code += self.__reset_element_tensor(tensors[0], format)

        # Loop tensors to generate quadrature loops
        for i in range(num_tensors):
            # Loop all quadrature points
            element_code += [indent(format_comment("Loop quadrature points (tensor %d)" %(i,)), self.indent_size)]
            element_code += self.__begin_loop("ip", len(tensors[i].quadrature.weights), format)

            # Generate code to evaluate the Jacobian at every quadrature point
            # FIXME: get dim and num_dofs more clever
            # FIXME: Switched off until a nodemap is available (assuming affine map)
#            element_code += self.__generate_jacobian(len(tensors[i].Derivatives), len(tensors[i].Psis[0][0][:,0]), i, #format)
            print "generating code for tensor :", i
            # Generate the element tensor
            element_code += self.__element_tensor(tensors[i], i, format)

            # End the quadrature loop
            element_code += self.__end_loop(format)

            if i + 1 < len(tensors):
                element_code += [""]

        return tabulate_code + element_code

    def __tabulate_weights(self, weights, tensor_number, format):
        "Generate table of quadrature weights"

        code = []
        
        # Prefetch formats to speed up code generation
        format_floating_point = format["floating point"]

        # Get number of weights
        num_weights = len(weights)

        code += [format["comment"]("Array of quadrature weights (tensor %d)" %(tensor_number,) )]

        # Create variable name
#        name = "const static double "+ format["weights"](tensor_number, str(num_weights))
        name = format["table declaration"] + format["weights"](tensor_number, str(num_weights))

        # Generate entries
        value = format["block begin"]
        for i in range(num_weights):
#            value += "%f" %(weights[i],)
            value += format_floating_point(weights[i])
            if i == num_weights - 1:
                value += format["block end"]
            else:
                value += ", "

        code += [(name, value)] + [""]

        return code

    def __tabulate_derivatives(self, Derivatives, tensor_number, format):
        "Generate table of derivatives at quadrature points"

        code = []
#        print "\n quad codegenerator, Jacobian: ", Derivatives
        
        # Get number of directional derivatives
        dim = len(Derivatives)
#        print "\n quad codegenerator, dim: ", dim


        if dim == 2:
            directions = [(1,0),(0,1)]
            # Loop directions
            for d in range(len(directions)):
                code += [format["comment"]("Table of derivatives at quadrature points in x%d-direction" %(d,) )]
                code += [format["comment"]("Format: [quadrature points][dofs] (tensor %d)" %(tensor_number,) )]

#                print "\n quad codegenerator, Derivatives[(1,0)]: ", Derivatives[directions[d]]

                # Get derivatives
                derivatives = Derivatives[directions[d]]
                print "derivatives: ", derivatives
                print numpy.shape(derivatives)

                # Get number of dofs (rows)
                num_dofs = len(derivatives[:,0])

                # Get number of quadrature points (columns)
                num_quadrature_points = len(derivatives[0,:])
                print "dofs: ", num_dofs
                print "ip: ", num_quadrature_points


                # Create variable name, static double array[num_dofs][num_quadraturepoints]
                # Create variable name, static double array[num_quadraturepoints][num_dofs]
                # Note that this is the transpose of what is returned from FIAT
#                name = "static double "+ format["derivatives"](d, str(num_dofs), str(num_quadrature_points))
                name = format["table declaration"] + format["derivatives"]\
                (tensor_number, d, str(num_quadrature_points), str(num_dofs))

                # Generate arrays
                value = "\\\n" + format["block begin"]
#                for i in range(num_dofs):
                for j in range(num_quadrature_points):
                    if j == 0:
                        value += format["block begin"]
                    else:
                        value += " " + format["block begin"]
#                    for j in range(num_quadrature_points):
                    for i in range(num_dofs):
#                        value += "%f" %(derivatives[i,j],)
                        value += format["floating point"] %(derivatives[i,j],)
#                        if j == num_quadrature_points - 1:
                        if i == num_dofs - 1:
                            value += format["block end"]
                        else:
                            value += ", "
#                    if i == num_dofs - 1:
                    if j == num_quadrature_points - 1:
                        value += format["block end"]
                    else:
                        value += ",\n"

                code += [(name, value)] + [""]

        # FIXME: Implement 3D support
        else:
            RuntimeError("tabulation of derivatives not implemented for 3D")

        return code

    def __tabulate_psis(self, tensor, tensor_number, format):
        "Tabulate values of basis functions and their derivatives at quadrature points"

        code = []

        # Get list of psis
        psis = tensor.Psis

        # Prefetch formats to speed up code generation
        format_psis = format["psis"]
        format_floating_point = format["floating point"]
        format_block_begin = format["block begin"]
        format_block_end = format["block end"]

#        print "\n quad codegenerator, tabulate_psis: ", psis

        # Get number of psis to tabulate
        num_psis = len(psis)
#        print "\n quad codegenerator, num_psis: ", num_psis

        code += [format["comment"]("Values of shapefunctions and their derivatives at quadrature points")]
        code += [format["comment"]("Format: [quadrature points][dofs] (tensor %d)" % (tensor_number,))]

        # Loop psis
        for psi_number in range(num_psis):

            # Get values of psi
            psi = psis[psi_number][0]

#            psi = psis[d]
#            print "psi: ", psi
#            print "psis[d][1]: ", psis[d][1]
#            print "psis[d][1]: ", psis[d][1]
#            print "numpy.shape(psi): ", numpy.shape(psi)
#            print "tensor.a: ", tensor.a
#            rank = tensor.a.rank

            # Get list of indices
            indices = psis[psi_number][1]

#            print "indices: ", indices
#            print "tuple(indices[0]): ", tuple(indices[0])
#            print "psis[tuple(indices[0])]: ", psi[tuple(indices[0])]
#            aindices = [[a] for a in range(rank)]
#            print "psis[d][2][0]: ", psis[d][2][0]

            # Get number of dofs, firste entry in shape list
            num_dofs = psis[psi_number][2][0]

#            print "len(tensor.quadrature.weights): ", len(tensor.quadrature.weights)

            # Get number of quadrature points
            num_quadrature_points = len(tensor.quadrature.weights)

            # Get index information for name generation
#            index_type = 0
#            index_index = 0
#            dims = []
#            for i in range(len(indices)):
#                index = indices[i]
#                print "index: ", index
#                print "index.type: ", index.type
#                print "index.index: ", index.index
#                if i == 0:
#                    index_type = index.type
#                    index_index = index.index
#                else:
#                    dims += [psis[psi_number][2][i]]
#                    print "dims: ", dims

            # Generate multi indices from the dimension(s) of the secondary indices of the current psi
#            aindices = MultiIndex(dims).indices

#            print "MultiIndex(dims).indices: ", aindices

            # Generate list of names
#            names = []
#            for multiindex in aindices:
#                print "multi index: ", multiindex
#                names += ["const static double " + format_psis([tensor_number, index_type, \
#                                index_index, multiindex, num_quadrature_points, num_dofs])]

            (names, aindices) = self.__generate_psi_declaration(tensor_number, indices,\
                                psis[psi_number][2], num_quadrature_points, num_dofs, format)

#            print "names: ", names

#            print "aindices: ", aindices
#            print "len(aindices): ", len(aindices)
#            print "tuple(aindices[0]): ", tuple(aindices[0])

            # Loop secondary indices
            for aindex_number in range(len(aindices)):
#            for a in range(rank) or [[]]:
#                print "a: ", a
#                print "tuple(a): ", tuple(a)
#                print "list(a): ", list(a)

                # Get name
                name = names[aindex_number]

                # Get list of secondary indices
                aindex = aindices[aindex_number]

                # Get values from psi tensor, should have format values[dofs][quad_points]
                values = psi[tuple(aindex)]

#                print "values: ", values
#            print "len(numpy.shape(psi)): ", len(numpy.shape(psi))
                # Get values of dofs at quadrature points
#                values = psi[a]
#                values = psi[tuple(a)]
#                print "values: ", values
#                shape = numpy.shape(values)
#                print "shape: ", shape

#                print "num_dofs: ", num_dofs

#                print "num_quadrature_points: ", num_quadrature_points

                # Create variable name, const static double name[num_quadraturepoints][num_dofs]
#                name = "const static double " + format_psis(tensor_number, d, a, \
#                                                str(num_quadrature_points), str(num_dofs))

                # Generate array of values
                value = "\\\n" + format_block_begin
                for i in range(num_quadrature_points):
                    if i == 0:
                        value += format_block_begin
                    else:
                        value += " " + format_block_begin
                    for j in range(num_dofs):
                        value += format_floating_point(values[j,i])
                        if j == num_dofs - 1:
                            value += format_block_end
                        else:
                            value += ", "
                    if i == num_quadrature_points - 1:
                        value += format_block_end
                    else:
                        value += ",\n"

                code += [(name, value)] + [""]
#        print "code: ", code
        return code

    def __generate_jacobian(self, dim, num_dofs, tensor_number, format):
        "Generate Jacobian at integration points"

        # Prefetch formats to speed up code generation
        format_comment        = format["comment"]
        format_transform      = format["transform"]
        format_determinant    = format["determinant"]
        format_derivatives    = format["derivatives"]
        format_coordinates    = format["coordinates"]
        format_add_equal      = format["add equal"]
        format_add            = format["add"]
        format_subtract       = format["subtract"]
        format_multiply       = format["multiply"]
        format_floating_point = format["floating point"]
        format_float = format["float declaration"]

        code = []

        # Declare variables
        code += [indent(format_comment("Declare and initialize variables for Jacobian and the determinant (tensor %d)" % (tensor_number,)), self.indent_size)]

        # FIXME: hardcoded variables
        code += [(indent(format_float + format_transform(0,0,None), self.indent_size), format_floating_point(0.0))]
        code += [(indent(format_float + format_transform(0,1,None), self.indent_size), format_floating_point(0.0))]
        code += [(indent(format_float + format_transform(1,0,None), self.indent_size), format_floating_point(0.0))]
        code += [(indent(format_float + format_transform(1,1,None), self.indent_size), format_floating_point(0.0))]
        code += [(indent(format_float + format_determinant, self.indent_size), format_floating_point(0.0))]
        code += [""]

        # Create loop over dofs
        code += [indent(format_comment("Jacobian, loop dofs (tensor %d)" % (tensor_number,)), self.indent_size)]
        code += self.__begin_loop("j", num_dofs, format)

        # Compute Jacobian values and determinant
#        code += [format_comment("Compute values of the Jacobian and determinant, loop dofs", self.indent_size)]
        if dim == 2:
            for i in range(dim):
                for j in range(dim):
                    # FIXME: hardcoded choose_map[], (3rd argument)
                    name = indent(format_transform(i,j, None), self.indent_size)
#                    print "name: ", name
#                    value = "dNdx%d[%s][%s]" % (i,"j","i")
                    value = format_multiply([format_derivatives(tensor_number, i,"ip","j"), format_coordinates("j",str(j))])
#                    print "value: ", value

#                    print "format_add_equal: ", format_add_equal(1,1)
#                    print "format_add_equal(name, value): ", format_add_equal(name,value)

                    code += [format_add_equal(name, value)]

#                    code += [(name, value)]

            code += [""]
            code += [indent(format_comment("Compute determinant of Jacobian"), self.indent_size)]

            # FIXME: determinant is hardcoded
            value = format_subtract([format_multiply([format_transform(0,0,None),format_transform(1,1,None)]), \
                                     format_multiply([format_transform(0,1,None),format_transform(1,0,None)])])

            code += [(indent(format_determinant, self.indent_size), value)]
            code += [""]

            code += [indent(format_comment("Take absolute value of determinant"), self.indent_size)]
            code += [(indent(format_determinant, self.indent_size), format["absolute value"](format_determinant))]

        else:
            RuntimeError("Jacobian for 3D not implemented yet!")

        # End node loop
        code += self.__end_loop(format)
        code += [""]

        return code

    def __element_tensor(self, tensor, tensor_number, format):
        "Generate loop over primary indices"

        code = []

        # Prefetch formats to speed up code generation
        format_element_tensor = format["element tensor"]
#        format_comment        = format["indent comment"]
        format_comment        = format["comment"]
        format_determinant    = format["determinant"]
        format_weights        = format["weights"]
        format_add_equal      = format["add equal"]
        format_add            = format["add"]
        format_subtract       = format["subtract"]
        format_multiply       = format["multiply"]
        format_floating_point = format["floating point"]
        format_psis           = format["psis"]

        # Get number of dofs, assuming all psis have the same number of dofs.
        num_dofs = tensor.Psis[0][2][0]

        # Loop first primary index
#        code += [format_comment("Loop first primary index", self.indent_size)]
        code += [indent(format_comment("Loop first primary index"), self.indent_size)]
        code += self.__begin_loop("i", num_dofs, format)

        # Loop second primary index
#        code += [format_comment("Loop second primary index", self.indent_size)]
        code += [indent(format_comment("Loop second primary index"), self.indent_size)]
        code += self.__begin_loop("j", num_dofs, format)

        scaling = [format_multiply([format_weights(tensor_number, "ip"), format_determinant]) + \
                  " \\\n%s" %("".join([" " for i in range(self.indent_size)]))]
#        print "scaling: ", scaling
#        code += [format_comment("Compute block entries (tensor %d)" % (tensor_number,), self.indent_size)]
        code += [indent(format_comment("Compute block entries (tensor %d)" % (tensor_number,)), self.indent_size)]

        iindices = tensor.i.indices
        aindices = tensor.a.indices

        indices = [psi[1] for psi in tensor.Psis]
#        print "indices: ", indices

        # Generate entry and name
        k = "i*%d+j" %(num_dofs,)
        name = indent(format["element tensor quad"](k), self.indent_size)

        value = []
        for a in aindices:
#            print "a"
            factor = self.__generate_factor(tensor, a, format)
            value += [format_multiply([self.__generate_psi_entry(tensor_number, a, psi_indices, format) for psi_indices in indices] + factor + scaling)]

        value = format_add(value)

        code += [format_add_equal(name, value)]
#        code += [(name, value)]

        # End loop second primary index
        code += self.__end_loop(format)

        # End loop first primary index
        code += self.__end_loop(format)

#        for i in iindices:
#            name = indent(format_element_tensor(i, k), self.indent_size)
#            value = []

#            for a in aindices:
#                factor = self.__generate_factor(tensor, a, format)
#                value += [format_multiply([self.__generate_psi_entry(tensor_number, i, a, ia, format) for ia in #indices] + factor + scaling)]


#            k += 1

        return code

    def __reset_element_tensor(self, tensor, format):
        "Reset the entries of the local element tensor"

        code = []

        # Prefetch formats to speed up code generation
        format_comment = format["comment"]

        # Get number of dofs, assuming all psis have the same number of dofs.
        num_dofs = tensor.Psis[0][2][0]

        # Loop first primary index
        code += [indent(format_comment("Reset values of the element tensor block, loop first primary index"), self.indent_size)]
        code += self.__begin_loop("i", num_dofs, format)

        # Loop second primary index
        code += [indent(format_comment("Loop second primary index"), self.indent_size)]
        code += self.__begin_loop("j", num_dofs, format)

        # Generate entry, name and value
        k = "i*%d+j" %(num_dofs,)
        name = indent(format["element tensor quad"](k), self.indent_size)
        value = format["floating point"](0.0)

        code += [indent(format_comment("Reset entry"), self.indent_size)]
        code += [(name, value)]

        # End loop second primary index
        code += self.__end_loop(format)

        # End loop first primary index
        code += self.__end_loop(format)

        code += [""]

        return code

    def __generate_psi_declaration(self, tensor_number, psi_indices, psi_shapes,\
                                         num_quadrature_points, num_dofs, format):

        psi_args = [tensor_number, 0, 0, [], num_quadrature_points, num_dofs]
        dims = []
        for i in range(len(psi_indices)):
            index = psi_indices[i]
            if i == 0:
                psi_args[1] = index.type
                psi_args[2] = index.index
            else:
                dims += [psi_shapes[i]]

        # Generate multi indices from the dimension(s) of the secondary indices of the current psi
        aindices = MultiIndex(dims).indices

        # Generate list of names
        names = []
        for multiindex in aindices:
            psi_args[3] = multiindex
            names += [format["table declaration"] + format["psis"](psi_args)]
 
        return (names, aindices)

#    def __generate_psi_entry(self, tensor_number, primary_indices, secondary_indices, psi_indices, format):
    def __generate_psi_entry(self, tensor_number, secondary_indices, psi_indices, format):

        psi_args = [tensor_number, 0, 0, [], "ip", 0]

        primary_indices = ["i", "j"]
#        print "psi_indices: ", psi_indices

        for i in range(len(psi_indices)):
            index = psi_indices[i]
            if i == 0:
                psi_args[1] = index.type
#                print "index.type: ", index.type
                psi_args[2] = index.index
#                print "index.index: ", index.index

                psi_args[5] = index(primary_indices, secondary_indices, [], [])
#                print "index(primary_indices, secondary_indices, [], []) :", index(primary_indices, secondary_indices, [], [])

            else:
                psi_args[3] += [index(primary_indices, secondary_indices, [], [])]
#                print "psi_args[3]: ", psi_args[3]
#        print "psi_args: ", psi_args

        return format["psis"](psi_args)

    def __generate_factor(self, tensor, multiindex, format):
        "Generate code for the value of entry a of geometry tensor G"

        # Compute product of factors outside sum
        factors = []
        for c in tensor.constants:
            if c.inverted:
                factors += ["(1.0/" + format["constant"](c.number.index) + ")"]
            else:
                factors += [format["constant"](c.number.index)]
        for c in tensor.coefficients:
            if not c.index.type == Index.AUXILIARY_G:
#                print "coefficient: ", c
#                print "multiindex: ", multiindex
#                print "c.n1.index: ", c.n1.index
#                print "c.index: ", c.index
#                print "eval multiindex: ", c.index([], multiindex, [], [])
                coefficient = format["coefficient"](c.n1.index, c.index([], multiindex, [], []))
                factors += [coefficient]
        for t in tensor.transforms:
#            print "t: ", t
            if not (t.index0.type == Index.AUXILIARY_G or  t.index1.type == Index.AUXILIARY_G):
#                print "transform: ", t
#                print "transform power: ", t.power
#                print "multiindex: ", multiindex
#                print "i0: ", t.index0
#                print "i1: ", t.index1
#                print "eval i0: ", t.index0([], multiindex, [], [])
#                print "eval i1: ", t.index1([], multiindex, [], [])
                factors += [format["inverse transform"](t.index0([], multiindex, [], []), \
                                                        t.index1([], multiindex, [], []), \
                                                        t.restriction),]
        return factors

#        monomial = format["multiply"](factors)
#        if monomial: f0 = [monomial]
#        else: f0 = []
    
        # Compute sum of monomials inside sum
#        terms = []
#        for b in G.b.indices:
#            factors = []
#            for c in G.coefficients:
#                if c.index.type == Index.AUXILIARY_G:
#                    coefficient = format["coefficient"](c.n1.index, c.index([], a, [], b))
#                    factors += [coefficient]
#            for t in G.transforms:
#                if t.index0.type == Index.AUXILIARY_G or t.index1.type == Index.AUXILIARY_G:
#                    factors += [format["inverse transform"](t.index0([], a, [], b), \
#                                                            t.index1([], a, [], b), \
#                                                            t.restriction)]
#            terms += [format["multiply"](factors)]
#        sum = format["add"](terms)
#        if sum: sum = format["grouping"](sum)
#        if sum: f1 = [sum]
#        else: f1 = []
    
        # Compute product of all factors
#        return format["multiply"]([f for f in [format["scale factor"]] + f0 + f1])

    def __begin_loop(self, variable, num_loops, format):
        "Generate the beginning of a loop"

        code = [indent(format["loop"](variable, variable, num_loops, variable), self.indent_size)]
        code += [indent(format["block begin"], self.indent_size)]

        # Increment indentation size
        self.indent_size += self.indent_increment

        return code

    def __end_loop(self, format):
        "Ends a loop"

        # Decrement indentation size
        self.indent_size -= self.indent_increment

        return [indent(format["block end"], self.indent_size)]



