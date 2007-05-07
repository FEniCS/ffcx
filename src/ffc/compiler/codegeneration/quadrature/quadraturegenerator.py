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

# FFC code generation modules
from ffc.compiler.codegeneration.common.codegenerator import *
from ffc.compiler.codegeneration.common.utils import *

from ffc.compiler.representation.quadrature.multiindex import *

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

    def __generate_element_tensor(self, tensors, Indent, format):
        "Construct quadrature code for element tensors"

        # Initialize code segments
        tabulate_code = []
        element_code = []

        # Prefetch formats to speed up code generation
        format_comment      = format["comment"]
        format_loop         = format["loop"]
        format_ip           = format["loop integration points"]
        format_block_begin  = format["block begin"]
        format_block_end    = format["block end"]

        # Number of tensors to evaluate
        num_tensors = len(tensors)

        # Loop tensors to generate tables
        for i in range(num_tensors):

            # Tabulate variables:

            # Tabulate the derivatives of map, used to generate the Jacobian
            # FIXME: Switched off until a nodemap is available (Jacobian), how to specify
            # on which element the map is defined???
#            tabulate_code += self.__generate_map_coordinates(tensors[i].map_element, Indent, format)
#            tabulate_code += self.__tabulate_derivatives(tensors[i].map_derivatives, i, Indent, format)

            # Tabulate the quadrature weights
            tabulate_code += self.__tabulate_weights(tensors[i].quadrature.weights, i, Indent, format)

            # Tabulate values of basis functions and their derivatives at quadrature points
            tabulate_code += self.__tabulate_psis(tensors[i], i, Indent, format)

        # Reset values of the element tensor (assuming same dimensions for all tensors)
        tabulate_code += self.__reset_element_tensor(tensors[0], Indent, format)

        # Loop tensors to generate quadrature loops
        for i in range(num_tensors):
            # Loop all quadrature points
            element_code += [Indent.indent(format_comment("Loop quadrature points (tensor/monomial term %d)" %(i,)))]
            element_code += [Indent.indent(format_loop(format_ip, 0, len(tensors[i].quadrature.weights)))]
            element_code += [Indent.indent(format_block_begin)]

            # Increase indentation
            Indent.increase()

            # Generate code to evaluate the Jacobian at every quadrature point
            # FIXME: Switched off until a nodemap is available (assuming affine map)
#            element_code += self.__generate_jacobian(tensors[i].map_element.cell_shape(),\
#                                                     tensors[i].map_element.space_dimension(), i, Indent, format)
            print "generating code for tensor :", i

            # Generate the element tensor
            element_code += self.__element_tensor(tensors[i], i, Indent, format)

            # Decrease indentation
            Indent.decrease()
            # End the quadrature loop
            element_code += [Indent.indent(format_block_end)]

            if i + 1 < len(tensors):
                element_code += [""]

        return tabulate_code + element_code

    def __tabulate_derivatives(self, map_derivatives, tensor_number, Indent, format):
        "Generate table of the first derivatives at quadrature points, used for Jacobian."

        code = []

        # Prefetch formats to speed up code generation
        format_comment      = format["comment"]
        format_table        = format["table declaration"]
        format_derivatives  = format["derivatives"]
        
        # Get number of directional derivatives
        dim = len(map_derivatives)

        if dim == 2:
            # Put this into a dictionary
            directions = [(1,0),(0,1)]
            # Loop directions
            for d in range(len(directions)):
                code += [Indent.indent(format_comment\
                        ("Table of derivatives of map at quadrature points in x%d-direction" %(d,) ))]
                code += [Indent.indent(format_comment\
                        ("Format: [quadrature points][dofs] (tensor %d)" %(tensor_number,) ))]

                # Get derivatives
                derivatives = map_derivatives[directions[d]]

                # Get number of dofs (rows)
                num_dofs = len(derivatives[:,0])

                # Get number of quadrature points (columns)
                num_quadrature_points = len(derivatives[0,:])

                # Create variable name, static double array[num_dofs][num_quadraturepoints]
                # Create variable name, static double array[num_quadraturepoints][num_dofs]
                # Note that this is the transpose of what is returned from FIAT
                name = format_table + format_derivatives\
                       (tensor_number, d, str(num_quadrature_points), str(num_dofs))

                # Generate arrays
                value = tabulate_matrix(numpy.transpose(derivatives), format)

                code += [(Indent.indent(name), value)] + [""]

        # FIXME: Implement 3D support
        else:
            RuntimeError("tabulation of derivatives not implemented for 3D")

        return code

    def __tabulate_weights(self, weights, tensor_number, Indent, format):
        "Generate table of quadrature weights"

        code = []
        
        # Prefetch formats to speed up code generation
        format_floating_point = format["floating point"]

        # Get number of weights
        num_weights = len(weights)

        code += [Indent.indent(format["comment"]("Array of quadrature weights (tensor/monomial term %d)" %(tensor_number,) ))]

        # Create variable name
        name = format["table declaration"] + format["weights"](tensor_number, str(num_weights))
        value = format["block"](format["separator"].join([format_floating_point(weights[i])\
                 for i in range(num_weights)]))

        code += [(Indent.indent(name), value)]

        return code + [""]



    def __tabulate_psis(self, tensor, tensor_number, Indent, format):
        "Tabulate values of basis functions and their derivatives at quadrature points"

        code = []

        # Prefetch formats to speed up code generation
        format_comment = format["comment"]

        format_psis = format["psis"]
        format_floating_point = format["floating point"]
        format_block_begin = format["block begin"]
        format_block_end = format["block end"]

        # Get list of psis
        psis = tensor.Psis

#        print "\n quad codegenerator, tabulate_psis: ", psis
#        print "\n quad codegenerator, tabulate_psis[0]: ", psis[0]
#        print "\n quad codegenerator, tabulate_psis[0][0]: ", psis[0][0]
        print "\n quad codegenerator, shape(tabulate_psis): ", numpy.shape(psis[0][0])

        # Get number of psis to tabulate
        num_psis = len(psis)
        print "\n quad codegenerator, num_psis: ", num_psis

        code += [Indent.indent(format_comment\
                ("Values of shapefunctions and their derivatives at quadrature points"))]
        code += [Indent.indent(format_comment\
                ("Format: [quadrature points][dofs] (tensor/monomial term %d)" % (tensor_number,)))]

        # Loop psis
        for psi_number in range(num_psis):

            # Get psi
            psi = psis[psi_number]

            # Get values of psi
            values = psi[0]

            # Get list of indices
            indices = psi[1]
            print "indices: ", indices
            # Get number of dofs, firste entry in shape list
            num_dofs = psi[2][0]

#            print "len(tensor.quadrature.weights): ", len(tensor.quadrature.weights)

            # Get number of quadrature points
            num_quadrature_points = len(tensor.quadrature.weights)

#            multiindices = psi[3]
#            print "multiindices: ", multiindices
#            for mul in multiindices:
#                print "mul: ", mul

            aindices = tensor.a.indices
            bindices = tensor.b.indices
#            aindices = psi[3][0]
#            bindices = psi[3][1]

            names = []
            multi_indices = []
            print "aindices: ", aindices
            print "bindices: ", bindices

            for a in aindices:
                for b in bindices:

                    (name, multi_index) = self.__generate_psi_declaration(tensor_number, indices,\
                                a, b, num_quadrature_points, num_dofs, format)

                    names += [name]
                    multi_indices += [multi_index]

            # Remove redundant names and entries in the psi table
            unique_names = []
            unique_multi_indices = []
            for i in range(len(names)):
                name = names[i]
                if not name in unique_names:
                    unique_names += [name]
                    unique_multi_indices += [multi_indices[i]]

            for i in range(len(unique_names)):
                # Get values from psi tensor, should have format values[dofs][quad_points]
                vals = values[tuple(unique_multi_indices[i])]

                # Generate array of values (FIAT returns [dof, ip] transpose to [ip, dof])
                value = tabulate_matrix(numpy.transpose(vals), format)

                code += [(Indent.indent(unique_names[i]), Indent.indent(value))] + [""]

                    # Get values from psi tensor, should have format values[dofs][quad_points]
#                    vals = values[tuple(multi_index)]

                    # Generate array of values (FIAT returns [dof, ip] transpose to [ip, dof])
#                    value = tabulate_matrix(numpy.transpose(vals), format)

#                    code += [(Indent.indent(name), Indent.indent(value))] + [""]


#            (names, aindices) = self.__generate_psi_declaration(tensor_number, indices,\
#                                psis[psi_number][2], num_quadrature_points, num_dofs, format)

#            print "aindices: ", aindices
            # Loop secondary indices
#            for aindex_number in range(len(aindices)):

                # Get name
#                name = names[aindex_number]

                # Get list of secondary indices
#                aindex = aindices[aindex_number]
#                print "aindex: ", aindex
#                print "tuple(aindex): ", tuple(aindex)
                # Get values from psi tensor, should have format values[dofs][quad_points]
#                values = psi[tuple(aindex)]

                # Generate array of values (FIAT returns [dof, ip] transpose to [ip, dof])
#                value = tabulate_matrix(numpy.transpose(values), format)

#                code += [(Indent.indent(name), Indent.indent(value))] + [""]

        return code

    def __reset_element_tensor(self, tensor, Indent, format):
        "Reset the entries of the local element tensor"

        code = []

        # Prefetch formats to speed up code generation
        format_comment  = format["comment"]
        format_i         = format["first free index"]
        format_j         = format["second free index"]

        # Get primary indices (MultiIndex)
        primary_indices = tensor.i

        # Generate value
        value = format["floating point"](0.0)

        # FIXME: quadrature only supports Linear and Bilinear forms
        if (primary_indices.rank == 1):
            code += [Indent.indent(format_comment\
                    ("Reset values of the element tensor block"))]

            # Generate name
            name =  format["element tensor quad"] + format["array access"](format_i)

            # Create boundaries for loop
            boundaries = [0, primary_indices.dims[0]]
            code += self.__generate_loop(name, value, boundaries, Indent, format)

        elif (primary_indices.rank == 2):
            code += [Indent.indent(format_comment\
                    ("Reset values of the element tensor block"))]

            # Generate entry and name
            entry = format["add"]([format["multiply"]([format_i, "%d" %primary_indices.dims[1]]), format_j])
            name =  format["element tensor quad"] + format["array access"](entry)

            # Create boundaries for loop
            boundaries = [0, primary_indices.dims[0], 0, primary_indices.dims[1]]
            code += self.__generate_loop(name, value, boundaries, Indent, format)
        else:
            raise RuntimeError, "Quadrature only supports Linear and Bilinear forms"

        return code + [""]

    def __generate_jacobian(self, dim, num_dofs, tensor_number, Indent, format):
        "Generate Jacobian at integration points"

        # Prefetch formats to speed up code generation
        format_comment        = format["comment"]
        format_float          = format["float declaration"]
        format_transform      = format["transform"]
        format_floating_point = format["floating point"]
        format_determinant    = format["determinant"]
        format_derivatives    = format["derivatives"]
        format_coordinates    = format["argument coordinates"]
        format_matrix_access  = format["matrix access"]
        format_multiply       = format["multiply"]
        format_add_equal      = format["add equal"]
        format_subtract       = format["subtract"]
        format_abs            = format["absolute value"]
        format_ip             = format["loop integration points"]
        format_i              = format["first free index"]

        code = []

        # Declare variables
        code += [Indent.indent(format_comment\
                ("Declare and initialize variables for Jacobian and the determinant (tensor %d)" % (tensor_number,)))]

        # Build indices from dimension
        indices = build_indices([dim, dim])

        # FIXME: hardcoded variables (restriction)
        for index in indices:
            code += [(Indent.indent(format_float + format_transform(-1,index[0],index[1],None)),\
                                    format_floating_point(0.0))]
        code += [(Indent.indent(format_float + format_determinant), format_floating_point(0.0))]
        code += [""]

        # Create loop over dofs
        code += [Indent.indent(format_comment("Jacobian, loop dofs (tensor %d)" % (tensor_number,)))]

        code += [Indent.indent(format["loop"](format_i, 0, num_dofs))]
        code += [Indent.indent(format["block begin"])]

        # Increase indentation
        Indent.increase()

        # FIXME: Need to get coordinates differently
        for i in range(dim):
            for j in range(dim):
                # FIXME: hardcoded choose_map[], (4th argument)
                name = format_transform(-1, i, j, None)
                value = format_multiply([format_derivatives(tensor_number, i,format_ip, format_i),\
                                         format_coordinates + format_matrix_access(format_i, str(j))])

                code += [format_add_equal(Indent.indent(name), value)]

        # End node loop
        # Decrease indentation
        Indent.decrease()
        code += [Indent.indent(format["block end"])]

        # Compute Jacobian values and determinant
        if dim == 2:

            code += [""]
            code += [Indent.indent(format_comment("Compute determinant of Jacobian"))]

            # FIXME: determinant is hardcoded
            value = format_subtract([format_multiply([format_transform(-1,0,0,None),format_transform(-1,1,1,None)]), \
                                     format_multiply([format_transform(-1,0,1,None),format_transform(-1,1,0,None)])])

            code += [(Indent.indent(format_determinant), value)]
            code += [""]

            code += [Indent.indent(format_comment("Take absolute value of determinant"))]
            code += [(Indent.indent(format_determinant), format_abs(format_determinant))]

        else:
            RuntimeError("Jacobian for 3D not implemented yet!")

        code += [""]

        return code

    def __element_tensor(self, tensor, tensor_number, Indent, format):
        "Generate loop over primary indices"

        code = []

        # Prefetch formats to speed up code generation
        format_element_tensor = format["element tensor"]
        format_comment        = format["comment"]
        format_determinant    = format["determinant"]
        format_weights        = format["weights"]
        format_add_equal      = format["add equal"]
        format_add            = format["add"]
        format_subtract       = format["subtract"]
        format_multiply       = format["multiply"]
        format_floating_point = format["floating point"]
        format_psis           = format["psis"]
        format_ip             = format["loop integration points"]
        format_i              = format["first free index"]
        format_j              = format["second free index"]

        # Get primary indices (MultiIndex)
        primary_indices = tensor.i

        # Get Psi indices, list of primary and secondary indices e.g. [[i0, a0], [i1, a1]]
        indices = [psi[1] for psi in tensor.Psis]

        # Get indices (combinations) from the primary and secondary MultiIndex
        iindices = tensor.i.indices
        aindices = tensor.a.indices
        bindices = tensor.b.indices
        print "aindices: ", aindices
        print "bindices: ", bindices
        # Treating secondary and auxilliary indices the same way
#        abindices = []
#        if (tensor.a.indices != [[]]):
#            print "tensor.a.indices: ", tensor.a.indices
#            abindices +=  tensor.a.indices
#        elif (tensor.b.indices != [[]]):
#            print "tensor.b.indices: ", tensor.b.indices
#            abindices +=  tensor.b.indices
#        else:
#            abindices = [[]]
#        print "abindices: ", abindices

        # Compute scaling
        weight = [format_weights(tensor_number, format_ip)]
        # Generate value
        values = []
        for a in aindices:
            for b in bindices:
#            print "a: ", a
                factor = self.__generate_factor(tensor, a, b, format)

                values += [format_multiply([self.__generate_psi_entry(tensor_number, a, b, psi_indices, format)\
                       for psi_indices in indices] + weight + factor) + format["new line"]]

        value = format_add(values)

        # FIXME: quadrature only supports Linear and Bilinear forms
        if (primary_indices.rank == 1):

            # Generate name
            name =  format["element tensor quad"] + format["array access"](format_i)
            code += [Indent.indent(format_comment("Compute block entries (tensor/monomial term %d)" % (tensor_number,)))]

            # Create boundaries for loop
            boundaries = [0, primary_indices.dims[0]]
            code += self.__generate_loop(name, value, boundaries, Indent, format, format_add_equal)

        elif (primary_indices.rank == 2):

            # Generate entry and name
            entry = format["add"]([format["multiply"]([format_i, "%d" %primary_indices.dims[1]]), format_j])
            name =  format["element tensor quad"] + format["array access"](entry)

            code += [Indent.indent(format_comment("Compute block entries (tensor/monomial term %d)" % (tensor_number,)))]

            # Create boundaries for loop
            boundaries = [0, primary_indices.dims[0], 0, primary_indices.dims[1]]
            code += self.__generate_loop(name, value, boundaries, Indent, format, format_add_equal)
        else:
            raise RuntimeError, "Quadrature only supports Linear and Bilinear forms"

        return code

#    def __generate_psi_declaration(self, tensor_number, psi_indices, psi_shapes,\
#                                         num_quadrature_points, num_dofs, format):
    def __generate_psi_declaration(self, tensor_number, psi_indices, aindices, bindices,\
                                         num_quadrature_points, num_dofs, format):

        # Prefetch formats to speed up code generation
        format_table            = format["table declaration"]
        format_psis             = format["psis"]
        format_matrix_access    = format["matrix access"]
        format_secondary_index  = format["secondary index"]

        # Should be in dictionary!!
        index_names = {0: lambda i: "fix_%s" %(i), 1: lambda i: "prim_%s" %(i), 2: lambda i: "sec_%s" %(i),\
                 3: lambda i: "aux_%s" %(i), 4: lambda i: "func_%s" %(i), 5: lambda i: "proj_%s" %(i),\
                 6: lambda i: "cons_%s" %(i), 7: lambda i: "aux0_%s" %(i), 8: lambda i: "auxG_%s" %(i)}

#        psi_args = [tensor_number, 0, 0, [], num_quadrature_points, num_dofs]
#        dims = []
        multi_index = []
        entries = []

        indices = ""
        for i in range(len(psi_indices)):
            index = psi_indices[i]
            if i == 0:
#                psi_args[1] = index.type
#                psi_args[2] = index.index
                indices += format_secondary_index(index_names[index.type](index.index))
            else:
#                dims += [psi_shapes[i]]
                indices += format_secondary_index(index_names[index.type](index.index))
                entries += [index([], aindices, bindices, [])]
                if (index.type == Index.SECONDARY or index.type == Index.AUXILIARY_0):
                    multi_index += [index([], aindices, bindices, [])]

        # Generate multi indices from the dimension(s) of the secondary indices of the current psi
#        aindices = MultiIndex(dims).indices


        # Generate list of names
#        names = []
#        for multiindex in aindices:
#            psi_args[3] = multiindex
#            multi = format_secondary_index("multi_" + "".join(str(mul) for mul in multiindex))

#            names += [format_table + format_psis + format_secondary_index("term_%d" %tensor_number)\
#                  + indices + multi + format_matrix_access(num_quadrature_points, num_dofs)]

        multi = format_secondary_index("multi_" + "".join(str(mul) for mul in multi_index))

        name = format_table + format_psis + format_secondary_index("term_%d" %tensor_number)\
                  + indices + multi + format_matrix_access(num_quadrature_points, num_dofs)

 
        return (name, entries)
#        return (names, aindices)

    def __generate_psi_entry(self, tensor_number, aindices, bindices, psi_indices, format):


        # Prefetch formats to speed up code generation
        format_psis             = format["psis"]
        format_matrix_access    = format["matrix access"]
        format_secondary_index  = format["secondary index"]

        # Should be in dictionary!!
        index_names = {0: lambda i: "fix_%s" %(i), 1: lambda i: "prim_%s" %(i), 2: lambda i: "sec_%s" %(i),\
                 3: lambda i: "aux_%s" %(i), 4: lambda i: "func_%s" %(i), 5: lambda i: "proj_%s" %(i),\
                 6: lambda i: "cons_%s" %(i), 7: lambda i: "aux0_%s" %(i), 8: lambda i: "auxG_%s" %(i)}

#        psi_args = [tensor_number, 0, 0, [], format["loop integration points"], 0]

        primary_indices = [format["first free index"], format["second free index"]]

        indices = ""
        multiindices = ""
        for i in range(len(psi_indices)):
            index = psi_indices[i]
            if i == 0:
#                psi_args[1] = index.type
#                print "index.type: ", index.type
#                psi_args[2] = index.index
#                print "index.index: ", index.index
                indices += format_secondary_index(index_names[index.type](index.index))

                dof_num = index(primary_indices, aindices, bindices, [])
#                psi_args[5] = index(primary_indices, secondary_indices, [], [])
#                print "index(primary_indices, secondary_indices, [], []) :", index(primary_indices, secondary_indices, [], [])

            else:
#                psi_args[3] += [index(primary_indices, aindices, bindices, [])]
                if (index.type == Index.SECONDARY or index.type == Index.AUXILIARY_0):
                    multiindices += str(index(primary_indices, aindices, bindices, []))

                indices += format_secondary_index(index_names[index.type](index.index))

            multi = format_secondary_index("multi_") + multiindices

#                print "psi_args[3]: ", psi_args[3]
#        print "psi_args: ", psi_args

#        entry = format_psis + format_matrix_access(format["loop integration points"], dof_num)
        entry = format_psis + format_secondary_index("term_%d" %tensor_number)\
                  + indices + multi + format_matrix_access(format["loop integration points"], dof_num)

        return entry

    def __generate_loop(self, name, value, boundaries, Indent, format, connect = None):
        "This function generates a loop over a vector or matrix."

        code = []

        # Prefetch formats to speed up code generation
        format_comment  = format["comment"]
        format_loop     = format["loop"]
        format_i        = format["first free index"]
        format_j        = format["second free index"]

        if (len(boundaries) == 2):
            # Loop index
            code += [Indent.indent(format_loop(format_i, boundaries[0], boundaries[1]))]
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
            code += [Indent.indent(format_loop(format_i, boundaries[0], boundaries[1]))]
            code += [Indent.indent(format["block begin"])]

            # Increase indentation
            Indent.increase()

            # Loop second primary index
            code += [Indent.indent(format_loop(format_j, boundaries[2], boundaries[3]))]

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

    def __generate_factor(self, tensor, a, b, format):
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
                coefficient = format["coefficient"](c.n1.index, c.index([], a, b, []))
                factors += [coefficient]
        for t in tensor.transforms:
#            print "t: ", t
            if not (t.index0.type == Index.AUXILIARY_G or  t.index1.type == Index.AUXILIARY_G):
#                print "transform: ", t
#                print "transform power: ", t.power
#                print "multiindex: ", multiindex
#                print "i0: ", t.index0
#                print "i1: ", t.index1
#                print "i0.type: ", t.index0.type
#                print "i1.type: ", t.index1.type
#                print "eval i0: ", t.index0([], multiindex, [], [])
#                print "eval i1: ", t.index1([], multiindex, [], [])
#                print "power: ", t.power
                factors += [format["transform"](t.power, t.index0([], a, b, []), \
                                                        t.index1([], a, b, []), \
                                                        t.restriction),]
        if tensor.determinant:
            d0 = format["power"](format["determinant"], tensor.determinant)
            d = format["multiply"]([format["scale factor"], d0])
        else:
            d = format["scale factor"]

        factors += [d]

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

    def __generate_map_coordinates(self, map_element, Indent, format):

        code = []

        # This function should not be needed when non-affine mappings are made available,
        # then the cell that is integrated should carry the information about dof coordinates.
        # Largely copied from dof_map.py __tabulate_coordinates()
 
        # Prefetch formats to speed up code generation
        format_coordinates          = format["argument coordinates"]
        format_element_coordinates  = format["element coordinates"]
        format_matrix_access        = format["matrix access"]
        format_block                = format["block"]
        format_separator            = format["separator"]

        # Get coordinates of the dofs (on FIAT reference element)
        points = map_element.dual_basis().pts

        # check if we get some points from fem.dofmap.py
        if (points):

            # Get the cell shape
            cell_shape = map_element.cell_shape()

#            code += [format["comment"]("Get cell vertices")]
#            code += [format["get cell vertices"]]

            code += [format["comment"]("DOF coordinates for mapping, assuming affine mapping!!")]

            name = format["const float declaration"] + format_coordinates + format_matrix_access(map_element.space_dimension(), cell_shape)

            # Create linear Lagrange element for the transformation
            element = FiniteElement("Lagrange", shape_to_string[cell_shape], 1)

            # Tabulate values of basisfunctions (on FIAT element)
            table = element.tabulate(0, points)

            # Get matrix of values of basisfunctions at points (dof, values at dofs on linear element)
            transformed_values = numpy.transpose(table[0][(0,)*cell_shape])

            # Get shape of matrix
            shape_val = numpy.shape(transformed_values)

            # Generate array of values
            value = format["new line"] + format["block begin"]
            rows = []

            for i in range(map_element.space_dimension()):
                cols = []
                for j in range(cell_shape):

                    values = [transformed_values[i][k] for k in range(shape_val[1])]
                    symbols = [format_element_coordinates(k,j) for k in range(shape_val[1])]
                    cols += [inner_product(values, symbols, format)]
                rows += [format_block(format_separator.join(cols))]

            value += format["block separator"].join(rows)
            value += format["block end"]
            code += [(name, value)]

        else:
            code += [format["exception"]("not implemented")]

        code += [""]

        return code
