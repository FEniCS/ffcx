"Code generator for quadrature representation"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2007-03-16 -- 2007-06-19"
__copyright__ = "Copyright (C) 2007 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Anders Logg 2007

# Python modules
import numpy
from sets import Set


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

# FFC format modules
from ffc.compiler.format.removeunused import *

class QuadratureGenerator(CodeGenerator):
    "Code generator for for tensor representation"

    def __init__(self):
        "Constructor"

        # Initialize common code generator
        CodeGenerator.__init__(self)
        self.optimise_level = 3
        self.save_tables = False
        self.unique_tables = True

    def generate_cell_integral(self, form_data, form_representation, sub_domain, format):
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

        # Generate element code + set of used geometry terms + set of used signs
        element_code, members_code, trans_set, signs_set = self.__generate_element_tensor\
                                                     (tensors, change_signs, None, None, Indent, format)

        # Remove unused declarations
        sign_code = self.__remove_unused(sign_code, signs_set, format)

        # Get Jacobian snippet
        jacobi_code = [format["generate jacobian"](form_data.cell_dimension, Integral.CELL)]

        # Remove unused declarations
        code = self.__remove_unused(jacobi_code, trans_set, format)

        # Add geometry tensor declarations and sign code
        code += sign_code

        # Add element code
        code += [""] + [format["comment"]("Compute element tensor")]
        code += element_code

        return {"tabulate_tensor": code, "members":members_code}

    def generate_exterior_facet_integral(self, form_data, form_representation, sub_domain, format):
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

        num_facets = len(tensors)
        cases = [None for i in range(num_facets)]
        trans_set = Set()
        for i in range(num_facets):
            case = [format_block_begin]

            # Assuming all tables have same dimensions for all facets (members_code)
            c, members_code, t_set, s_set = self.__generate_element_tensor(tensors[i], False, i, None, Indent, format)
            case += c
            case += [format_block_end]
            cases[i] = case
            trans_set = trans_set | t_set

        # Get Jacobian snippet
        jacobi_code = [format["generate jacobian"](form_data.cell_dimension, Integral.EXTERIOR_FACET)]

        # Remove unused declarations
        common = self.__remove_unused(jacobi_code, trans_set, format)

        # Add element code
        common += [""] + [format["comment"]("Compute element tensor for all facets")]

        return {"tabulate_tensor": (common, cases), "constructor":"// Do nothing", "members":members_code}
    
    def generate_interior_facet_integral(self, form_data, form_representation, sub_domain, format):
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

        num_facets = len(tensors)
        cases = [[None for j in range(num_facets)] for i in range(num_facets)]
        trans_set = Set()
        for i in range(num_facets):
            for j in range(num_facets):
                case = [format_block_begin]

                # Assuming all tables have same dimensions for all facet-facet combinations (members_code)
                c, members_code, t_set, s_set = self.__generate_element_tensor(tensors[i][j],\
                                                False, i, j, Indent, format)
                case += c
                case += [format_block_end]
                cases[i][j] = case
                trans_set = trans_set | t_set

        # Get Jacobian snippet
        jacobi_code = [format["generate jacobian"](form_data.cell_dimension, Integral.INTERIOR_FACET)]

        # Remove unused declarations
        common = self.__remove_unused(jacobi_code, trans_set, format)

        # Add element code
        common += [""] + [format["comment"]("Compute element tensor for all facets")]

        return {"tabulate_tensor": (common, cases), "constructor":"// Do nothing", "members":members_code}

    def __generate_element_tensor(self, tensors, sign_changes, facet0, facet1, Indent, format):
        "Construct quadrature code for element tensors"

        # Initialize code segments
        tabulate_code = []
        element_code = []
        trans_set = Set()
        signs_set = Set()

        # Prefetch formats to speed up code generation
        format_comment      = format["comment"]
        format_loop         = format["loop"]
        format_ip           = format["integration points"]
        format_block_begin  = format["block begin"]
        format_block_end    = format["block end"]

        # Group tensors after number of quadrature points, and dimension of primary indices
        # to reduce number of loops
        group_tensors = equal_loops(tensors)
        tables = None
        name_map = {}
        # Get dictionary of unique tables, and the name_map
        if self.unique_tables:
            name_map, tables = unique_tables(tensors, format)

        # Generate load_table.h if tables should be saved.
        members_code = ""
        if self.save_tables:
            members_code = generate_load_table(tensors)

        for i in range(len(tensors)):
            # Tabulate the quadrature weights
            tabulate_code += self.__tabulate_weights(tensors[i].quadrature.weights, i, Indent, format)

        if self.save_tables:
            # Save psi tables instead of tabulating
            tabulate_code += save_psis(tensors, facet0, facet1, Indent, format, tables)
        else:
            # Tabulate values of basis functions and their derivatives at quadrature points
            tabulate_code += self.__tabulate_psis(tensors, Indent, format, tables)

        # Reset values of the element tensor (assuming same dimensions for all tensors)
        tabulate_code += self.__reset_element_tensor(tensors[0], Indent, format)

        for points in group_tensors:
            # Loop all quadrature points
            # Create list of tensors for comment
            ts = [group_tensors[points][idims][i] for idims in group_tensors[points]\
                                                  for i in range(len(group_tensors[points][idims]))]
            element_code += [Indent.indent(format_comment\
            ("Loop quadrature points (tensor/monomial terms %s)" %(str(tuple(ts)))))]
            element_code += [Indent.indent(format_loop(format_ip, 0, points))]
            element_code += [Indent.indent(format_block_begin)]

            # Increase indentation
            Indent.increase()

            # Get dictionary of primary indices
            prim_dic = group_tensors[points]

            # Generate loop over primary indices
            for idims in prim_dic:
                tensor_numbers = prim_dic[idims]
                indices = [format["first free index"], format["second free index"]]
                # Create loop variables
                loop_vars = [[indices[i], 0, idims[i]] for i in range(len(idims))]

                # Generate loop
                element_code += generate_loop("", "", loop_vars, Indent, format, "")
                Indent.decrease()
                element_code += [Indent.indent(format_block_begin)]
                Indent.increase()

                # Generate element tensors for all tensors with the current number of quadrature points
                for i in tensor_numbers:
                    e_code, t_set, s_set = self.__element_tensor(tensors[i], i, sign_changes, Indent, format, name_map)
#                    element_code += self.__element_tensor(tensors[i], i, sign_changes, Indent, format, name_map)
                    element_code += e_code
                    trans_set = trans_set | t_set
                    signs_set = signs_set | s_set

                # End loop primary indices
                # Decrease indentation
                Indent.decrease()
                element_code += [Indent.indent(format_block_end)]
                # Decrease indentation
                for i in range(len(idims) - 1):
                    Indent.decrease()

            # End the quadrature loop
            # Decrease indentation
            Indent.decrease()
            element_code += [Indent.indent(format_block_end)]

            if i + 1 < len(tensors):
                element_code += [""]
        return (tabulate_code + element_code, members_code, trans_set, signs_set)

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

    def __tabulate_psis(self, tensors, Indent, format, tables):
        "Tabulate values of basis functions and their derivatives at quadrature points"

        code = []

        format_psis = format["psis"]
        format_floating_point = format["floating point"]
        format_block_begin = format["block begin"]
        format_block_end = format["block end"]

        if tables:
            for names in tables:
                # Get value and save as an array
                vals = tables[names]
                # Generate array of values (FIAT returns [dof, quad_points] transpose to [quad_points, dof])
                value = tabulate_matrix(vals, format)
                code += [(Indent.indent(names), Indent.indent(value))]# + [""]
        else:
            for tensor_number in range(len(tensors)):
                tensor = tensors[tensor_number]
                tables = get_names_tables(tensor, tensor_number, format)

                code += [Indent.indent(format["comment"]\
                        ("Values of shapefunctions and their derivatives at quadrature points"))]
                code += [Indent.indent(format["comment"]\
                        ("Format: [quadrature points][dofs] (tensor/monomial term %d)" % (tensor_number,)))]
                for names in tables:
                    # Get value and save as an array
                    vals = tables[names]
                    # Generate array of values (FIAT returns [dof, quad_points] transpose to [quad_points, dof])
                    value = tabulate_matrix(vals, format)
                    code += [(Indent.indent(names), Indent.indent(value))]# + [""]
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

        # FIXME: quadrature only support Functionals and Linear and Bilinear forms
        if (irank == 0):
            code += [Indent.indent(format["comment"]("Reset value"))]

            # Generate entry and name
            entry = "0"
            name =  format["element tensor quad"] + format["array access"](entry)
            code += [(Indent.indent(name),value)]

        elif (irank == 1):
            code += [Indent.indent(format["comment"]\
                    ("Reset values of the element tensor block"))]

            # Generate entry and name
            entry = format["first free index"]
            name =  format["element tensor quad"] + format["array access"](entry)

            # Create boundaries for loop
            boundaries = [0, macro_idims[0]]
            loop_vars = [[format["first free index"]] + boundaries]
            code += generate_loop(name, value, loop_vars, Indent, format)

        elif (irank == 2):
            code += [Indent.indent(format["comment"]\
                    ("Reset values of the element tensor block"))]

            # Generate entry and name
            entry = format["add"]([format["multiply"]([format["first free index"], "%d" %macro_idims[1]]),\
                                                       format["second free index"]])
            name =  format["element tensor quad"] + format["array access"](entry)

            # Create boundaries for loop
            boundaries = [[0, macro_idims[0]], [0, macro_idims[1]]]
            loop_vars = [[format["first free index"]] + boundaries[0],\
                         [format["second free index"]] + boundaries[1]]

            code += generate_loop(name, value, loop_vars, Indent, format)
        else:
            raise RuntimeError, "Quadrature only support Functionals and Linear and Bilinear forms"

        return code + [""]

    def __element_tensor(self, tensor, tensor_number, sign_changes, Indent, format, name_map):
        "Generate loop over primary indices"

        code = []
        s_set = Set()

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
            values, secondary_loop, t_set = values_level_0(indices, vindices, aindices, b0indices,\
                                            bgindices, tensor, tensor_number, weight, format, name_map)
        elif self.optimise_level == 1:
            values, secondary_loop, t_set = values_level_1(indices, vindices, aindices, b0indices,\
                                            bgindices, tensor, tensor_number, weight, format, name_map)
        elif self.optimise_level == 2:
            values, secondary_loop, t_set = values_level_2(indices, vindices, aindices, b0indices,\
                                            bgindices, tensor, tensor_number, weight, format, name_map)
        elif self.optimise_level == 3:
            values, secondary_loop, t_set = values_level_3(indices, vindices, aindices, b0indices,\
                                            bgindices, tensor, tensor_number, weight, format, name_map)
        else:
            raise RuntimeError, "Optimisation level not implemented!"

        value = format_add(values)

        # FIXME: quadrature only support Functionals and Linear and Bilinear forms
        if (irank == 0):

            if sign_changes:
                value, s_set = add_sign(value, 0, [], format)

            # Entry is zero because functional is a scalar value
            entry = "0"
            # Generate name
            name =  format["element tensor quad"] + format["array access"](entry)
            code += [Indent.indent(format["comment"]\
                    ("Compute value (tensor/monomial term %d)" % (tensor_number,)))]

            # Create boundaries for loop
            loop_vars = secondary_loop
            if secondary_loop:
                code += generate_loop(name, value, loop_vars, Indent, format, format["add equal"])
            else:
                code += [format["add equal"](Indent.indent(name), value)]

        elif (irank == 1):
            # Generate entry
            for i in range(irank):
                for v in monomial.basisfunctions:
                    if v.index.type == Index.PRIMARY and v.index.index == i:
                        if v.restriction == Restriction.MINUS:
                            entry = format_add([dic_indices[i], str(idims[0])])
                        else:
                            entry = dic_indices[i]
                        break

            if sign_changes:
                value, s_set = add_sign(value, 0, [format["first free index"]], format)

            # Generate name
            name =  format["element tensor quad"] + format["array access"](entry)
            code += [Indent.indent(format["comment"]\
                    ("Compute block entries (tensor/monomial term %d)" % (tensor_number,)))]

            # Create boundaries for loop
            loop_vars = secondary_loop
            if secondary_loop:
                code += generate_loop(name, value, loop_vars, Indent, format, format["add equal"])
            else:
                code += [format["add equal"](Indent.indent(name), value)]

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

            if sign_changes:
                value, s_set = add_sign(value, 0, [format["first free index"], format["second free index"]], format)

            entry[0] = format_multiply([entry[0], str(macro_idims[1])])
            name =  format["element tensor quad"] + format["array access"](format_add(entry))

            code += [Indent.indent(format["comment"]\
                    ("Compute block entries (tensor/monomial term %d)" % (tensor_number,)))]

            # Create boundaries for loop
            loop_vars = secondary_loop
            if secondary_loop:
                code += generate_loop(name, value, loop_vars, Indent, format, format["add equal"])
            else:
                code += [format["add equal"](Indent.indent(name), value)]

        else:
            raise RuntimeError, "Quadrature only support Functionals and Linear and Bilinear forms"

        return (code, t_set, s_set)

    def __remove_unused(self, code, set, format):

        if code:
            # Generate body of code, using the format
            lines = format["generate body"](code)

            # Generate auxiliary code line that uses all members of the set (to trick remove_unused)
            line_set = format["add equal"]("A", format["multiply"](set))
            lines += "\n" + line_set

            # Remove unused Jacobi declarations
            code = remove_unused(lines)

            # Delete auxiliary line
            code = code.replace("\n" + line_set, "")

            return [code]
        else:
            return code
