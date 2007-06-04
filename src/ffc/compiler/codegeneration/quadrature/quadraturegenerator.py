"Code generator for quadrature representation"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2007-03-16 -- 2007-06-01"
__copyright__ = "Copyright (C) 2007 Kristian B. Oelgaard"
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

# Utility functions for quadraturegenerator
from quadraturegenerator_utils import *

class QuadratureGenerator(CodeGenerator):
    "Code generator for for tensor representation"

    def __init__(self):
        "Constructor"

        # Initialize common code generator
        CodeGenerator.__init__(self)
        self.optimise_level = 2
        self.save_tables = True

    def generate_cell_integral(self, form_representation, sub_domain, format):
        """Generate dictionary of code for cell integral from the given
        form representation according to the given format"""

        code = []

        # Object to control the code indentation
        Indent = IndentControl()

        # Extract tensors
        tensors = form_representation.cell_tensor
        if len(tensors) == 0:
            return None

        # Generate code for sign changes
        (sign_code, change_signs) = generate_signs(tensors, format)
        code += sign_code

        # Generate code for element tensor(s)
        code += [Indent.indent(format["comment"]("Compute element tensor"))]
        code += self.__generate_element_tensor(tensors, change_signs, False, False, Indent, format)

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
            case += self.__generate_element_tensor(tensors[i], False, i, False, Indent, format)
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
                case += self.__generate_element_tensor(tensors[i][j], False, i, j, Indent, format)
                case += [format_block_end]
                cases[i][j] = case

        return {"tabulate_tensor": (common, cases)}

    def __generate_element_tensor(self, tensors, sign_changes, facet0, facet1, Indent, format):
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

        # Group tensors after number of quadrature points, reduce number of loops
        group_tensors = equal_num_quadrature_points(tensors)

        # Generate load_table.h if tables should be saved.
        if self.save_tables:
            generate_load_table(tensors)
        for points in group_tensors:
            # Loop tensors to generate tables
            for i in group_tensors[points]:
                # Tabulate the quadrature weights
                tabulate_code += self.__tabulate_weights(tensors[i].quadrature.weights, i, Indent, format)

                if self.save_tables:
                    # Save psi tables instead of tabulating
                    tabulate_code += save_psis(tensors[i], i, facet0, facet1, Indent, format)
                else:
                    # Tabulate values of basis functions and their derivatives at quadrature points
                    tabulate_code += self.__tabulate_psis(tensors[i], i, Indent, format)

        # Reset values of the element tensor (assuming same dimensions for all tensors)
        tabulate_code += self.__reset_element_tensor(tensors[0], Indent, format)

        for points in group_tensors:
            # Loop all quadrature points
            element_code += [Indent.indent(format_comment("Loop quadrature points (tensor/monomial term %d)" %(i,)))]
            element_code += [Indent.indent(format_loop(format_ip, 0, points))]
            element_code += [Indent.indent(format_block_begin)]

            # Increase indentation
            Indent.increase()

            # Generate element tensors for all tensors with the current number of quadrature points
            for i in group_tensors[points]:
                element_code += self.__element_tensor(tensors[i], i, sign_changes, Indent, format)

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

        code += [Indent.indent(format["comment"]\
                ("Array of quadrature weights (tensor/monomial term %d)" %(tensor_number,) ))]

        # Create variable name
        name = format["table declaration"] + format["weights"](tensor_number, str(len(weights)))
        value = format["block"](format["separator"].join([format_floating_point(w)\
                 for w in weights]))

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

                    (name, multi_index) = generate_psi_declaration(tensor_number, indices, vindex,\
                                a, b, num_quadrature_points, num_dofs, format)

                    names += [name]
                    multi_indices += [multi_index]

            # Remove redundant names and entries in the psi table
            names, multi_indices = extract_unique(names, multi_indices)

            # Loop names and tabulate psis
            for i in range(len(names)):
                # Get values from psi tensor, should have format values[dofs][quad_points]
                vals = values[tuple(multi_indices[i])]

                # Check if the values have the correct dimensions, otherwise quadrature is not correctly
                # implemented for the given form!!
                if numpy.shape(vals) != (num_dofs, num_quadrature_points):
                    raise RuntimeError, "Quadrature is not correctly implemented for the given form!"

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

        # Get monomial and compute macro dimensions in case of restricted basisfunctions
        monomial = tensor.monomial
        macro_idims = compute_macro_idims(monomial, idims, irank)

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
            code += generate_loop(name, value, boundaries, Indent, format)

        elif (irank == 2):
            code += [Indent.indent(format["comment"]\
                    ("Reset values of the element tensor block"))]

            # Generate entry and name
            entry = format["add"]([format["multiply"]([format["first free index"], "%d" %macro_idims[1]]),\
                                                       format["second free index"]])
            name =  format["element tensor quad"] + format["array access"](entry)

            # Create boundaries for loop
            boundaries = [0, macro_idims[0], 0, macro_idims[1]]
            code += generate_loop(name, value, boundaries, Indent, format)
        else:
            raise RuntimeError, "Quadrature only supports Linear and Bilinear forms"

        return code + [""]

    def __element_tensor(self, tensor, tensor_number, sign_changes, Indent, format):
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

        # Get monomial and compute macro dimensions in case of restricted basisfunctions
        monomial = tensor.monomial
        macro_idims = compute_macro_idims(monomial, idims, irank)

        # Get Psi indices, list of primary and secondary indices e.g. [[i0, a0], [i1, a1]]
        indices = [psi[1] for psi in tensor.Psis]
        vindices = [psi[2] for psi in tensor.Psis]

        # Get secondary and auxiliary multiindices
        aindices, b0indices, bgindices = tensor.a.indices, tensor.b0.indices, tensor.bg.indices

        # Compute scaling
        weight = [format["weights"](tensor_number, format["integration points"])]

        # Choose level of optimisation
        if self.optimise_level == 0:
            values = values_level_0(indices, vindices, aindices, b0indices, bgindices,\
                                    tensor, tensor_number, weight, format)
        elif self.optimise_level == 1:
            values = values_level_1(indices, vindices, aindices, b0indices, bgindices,\
                                    tensor, tensor_number, weight, format)
        elif self.optimise_level == 2:
            values = values_level_2(indices, vindices, aindices, b0indices, bgindices,\
                                    tensor, tensor_number, weight, format)
        else:
            raise RuntimeError, "Optimisation level not implemented!"

        value = format_add(values)

        if sign_changes:
            value = add_sign(value, 0, [format["first free index"], format["second free index"]], format)

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
            code += [Indent.indent(format["comment"]\
                    ("Compute block entries (tensor/monomial term %d)" % (tensor_number,)))]

            # Create boundaries for loop
            boundaries = [0, idims[0]]
            code += generate_loop(name, value, boundaries, Indent, format, format["add equal"])

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

            code += [Indent.indent(format["comment"]\
                    ("Compute block entries (tensor/monomial term %d)" % (tensor_number,)))]

            # Create boundaries for loop
            boundaries = [0, idims[0], 0, idims[1]]
            code += generate_loop(name, value, boundaries, Indent, format, format["add equal"])
        else:
            raise RuntimeError, "Quadrature only supports Linear and Bilinear forms"

        return code
