"Code generator for quadrature representation"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2007-03-16 -- 2008-02-05"
__copyright__ = "Copyright (C) 2007-2008 Kristian B. Oelgaard"
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

# Utility and optimisation functions for quadraturegenerator
from quadraturegenerator_utils import *
from quadraturegenerator_optimisation import *
import reduce_operations

# FFC format modules
from ffc.compiler.format.removeunused import *

class QuadratureGenerator(CodeGenerator):
    "Code generator for for tensor representation"

    def __init__(self):
        "Constructor"

        # Initialize common code generator
        CodeGenerator.__init__(self)
        self.optimise_level = 11
        self.save_tables = False
        self.unique_tables = True

    def generate_cell_integral(self, form_data, form_representation, sub_domain, format):
        """Generate dictionary of code for cell integral from the given
        form representation according to the given format"""

        code = []

        # Object to control the code indentation
        Indent = IndentControl()

        # Extract terms for sub domain
        tensors = [term for term in form_representation.cell_tensor if term.monomial.integral.sub_domain == sub_domain]
        if len(tensors) == 0:
            return None


        # Generate element code + set of used geometry terms
        element_code, members_code, trans_set = self.__generate_element_tensor\
                                                     (tensors, None, None, Indent, format)

        # Get Jacobian snippet
        jacobi_code = [format["generate jacobian"](form_data.cell_dimension, Integral.CELL)]

        # Remove unused declarations
        code = self.__remove_unused(jacobi_code, trans_set, format)

        # Add element code
        code += [""] + [format["comment"]("Compute element tensor (using quadrature representation, optimisation level %d)" %self.optimise_level)]
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

        # Extract terms for sub domain
        tensors = [[term for term in t if term.monomial.integral.sub_domain == sub_domain] for t in form_representation.exterior_facet_tensors]
        if len(tensors) == 0:
            return None

        num_facets = len(tensors)
        cases = [None for i in range(num_facets)]
        trans_set = Set()
        for i in range(num_facets):
            case = [format_block_begin]

            # Assuming all tables have same dimensions for all facets (members_code)
            c, members_code, t_set = self.__generate_element_tensor(tensors[i], i, None, Indent, format)
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

        # Extract terms for sub domain
        tensors = [[[term for term in t2 if term.monomial.integral.sub_domain == sub_domain] for t2 in t1] for t1 in form_representation.interior_facet_tensors]
        if len(tensors) == 0:
            return None

        num_facets = len(tensors)
        cases = [[None for j in range(num_facets)] for i in range(num_facets)]
        trans_set = Set()
        for i in range(num_facets):
            for j in range(num_facets):
                case = [format_block_begin]

                # Assuming all tables have same dimensions for all facet-facet combinations (members_code)
                c, members_code, t_set = self.__generate_element_tensor(tensors[i][j], i, j, Indent, format)
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

    def __generate_element_tensor(self, tensors, facet0, facet1, Indent, format):
        "Construct quadrature code for element tensors"

        # Initialize code segments
        tabulate_code = []

        # Prefetch formats to speed up code generation
        format_comment      = format["comment"]
        format_loop         = format["loop"]
        format_ip           = format["integration points"]
        format_block_begin  = format["block begin"]
        format_block_end    = format["block end"]

        # Group tensors after number of quadrature points, and dimension of primary indices
        # to reduce number of loops
        group_tensors = equal_loops(tensors)

        # Reset QuadratureElement indices if optimisation level is below 5
        if self.optimise_level <= 4:
            for t in tensors:
                t.qei = []

        psi_name_map = {}
        psi_tables = {}
        non_zero_columns = {}
        # Get dictionary of unique tables, and the name_map
        if self.unique_tables or self.optimise_level >= 6:
            psi_name_map, psi_tables, non_zero_columns = unique_psi_tables(tensors, self.optimise_level, format)

#        print "psi_tables: ", psi_tables

        # Generate load_table.h if tables should be saved.
        members_code = ""
        if self.save_tables:
            members_code = generate_load_table(tensors)

        weight_name_map = {}
        if self.optimise_level >= 6:
            weight_name_map = unique_weight_tables(tensors, format)

        # Tabulate the quadrature weights
        tabulate_code += self.__tabulate_weights(tensors, Indent, format, weight_name_map)
#        tabulate_code += self.__tabulate_weights(tensors[i].quadrature.weights, i, Indent, format)

        cols_name_map = {}
        if self.save_tables:
            # Save psi tables instead of tabulating
            tabulate_code += save_psis(tensors, facet0, facet1, Indent, format, psi_tables)
        else:
            # Tabulate values of basis functions and their derivatives at quadrature points
            t_code, cols_name_map = self.__tabulate_psis(tensors, Indent, format, psi_tables, non_zero_columns)
            tabulate_code += t_code
        # Reset values of the element tensor (assuming same dimensions for all tensors)
        tabulate_code += self.__reset_element_tensor(tensors[0], Indent, format)

        element_code = []
        trans_set = Set()
        # Detour to hybrid representation
        if self.optimise_level >= 6:
            e_code, t_set = self.__element_tensor2(tensors, group_tensors, Indent, format, psi_name_map, weight_name_map, non_zero_columns, cols_name_map)
            element_code += e_code
            trans_set = trans_set | t_set
        # Standard quadrature representation
        else:
            tensor_ops_count = 0
            for points in group_tensors:
                ip_code = []
                # Loop all quadrature points
                # Create list of tensors for comment
                ts = [group_tensors[points][idims][i] for idims in group_tensors[points]\
                                                  for i in range(len(group_tensors[points][idims]))]
                ip_code += [Indent.indent(format_comment\
                ("Loop quadrature points (tensor/monomial terms %s)" %(str(tuple(ts)))))]
                ip_code += [Indent.indent(format_loop(format_ip, 0, points))]
                ip_code += [Indent.indent(format_block_begin)]

                # Increase indentation
                Indent.increase()

                # Get dictionary of primary indices
                prim_dic = group_tensors[points]

                # Generate loop over primary indices
                ip_ops_count = 0
                for idims in prim_dic:
                    tensor_code = []
                    tensor_numbers = prim_dic[idims]
                    indices = [format["first free index"], format["second free index"]]
                    # Create loop variables
                    loop_vars = [[indices[i], 0, idims[i]] for i in range(len(idims))]

                    # Generate loop
                    tensor_code += generate_loop("", "", loop_vars, Indent, format, "")
                    Indent.decrease()
                    tensor_code += [Indent.indent(format_block_begin)]
                    Indent.increase()

                    # Generate element tensors for all tensors with the current number of quadrature points
                    for i in tensor_numbers:
                        e_code, t_set, ops_count = self.__element_tensor(tensors[i], i, Indent, format, psi_name_map)
                        # The operations are done for each of the primary loops
                        for idim in idims:
                            ops_count *= idim
                        tensor_code += e_code
                        trans_set = trans_set | t_set
                        ip_ops_count += ops_count
                    # End loop primary indices
                    # Decrease indentation
                    Indent.decrease()
                    tensor_code += [Indent.indent(format_block_end)]

                    ip_code += tensor_code
                    # Decrease indentation
                    for i in range(len(idims) - 1):
                        Indent.decrease()

                # End the quadrature loop
                # Decrease indentation
                Indent.decrease()
                ip_code += [Indent.indent(format_block_end)]

                if i + 1 < len(tensors):
                    ip_code += [""]

                ip_code = [Indent.indent(format_comment\
                ("Number of operations to compute element tensor for following IP loop = %d" %(ip_ops_count*points) ))] + ip_code
                element_code += ip_code
                # The operations are done for each IP
                tensor_ops_count += ip_ops_count*points

            element_code = [Indent.indent(format_comment\
            ("Total number of operations to compute element tensor for all IP loops = %d" %(tensor_ops_count) ))] + element_code

        return (tabulate_code + element_code, members_code, trans_set)


    def __tabulate_weights(self, tensors, Indent, format, name_map):
        "Generate table of quadrature weights"

        code = []
        
        # Prefetch formats to speed up code generation
        format_floating_point = format["floating point"]

        for tensor_number in range(len(tensors)):
            if name_map:
                # Only tabulate values for the tensor number that is not mapped
                if tensor_number in name_map:
                    weights = tensors[tensor_number].quadrature.weights
                    # Get the list of tensors that this weight is valid for
                    tw = [tensor_number] + name_map[tensor_number]
                    tw.sort()
                    code += [Indent.indent(format["comment"]\
                    ("Array of quadrature weights (tensor/monomial terms %s)" %(str(tuple(tw)))))]
                    # Create variable name
                    name = format["table declaration"] + format["weights"](tensor_number, str(len(weights)))
                    value = format["block"](format["separator"].join([format_floating_point(w)\
                         for w in weights]))
                    code += [(Indent.indent(name), value)]
            # Tabulate ALL weights
            else:
                weights = tensors[tensor_number].quadrature.weights

                code += [Indent.indent(format["comment"]\
                    ("Array of quadrature weights (tensor/monomial term %d)" %(tensor_number,) ))]

                # Create variable name
                name = format["table declaration"] + format["weights"](tensor_number, str(len(weights)))
                value = format["block"](format["separator"].join([format_floating_point(w)\
                         for w in weights]))

                code += [(Indent.indent(name), value)]

        return code + [""]

    def __tabulate_psis(self, tensors, Indent, format, tables, non_zero_cols):
        "Tabulate values of basis functions and their derivatives at quadrature points"

        code = []
        cols_name_map = {}

        format_psis = format["psis"]
        format_floating_point = format["floating point"]
        format_block_begin = format["block begin"]
        format_block_end = format["block end"]

#        print "non_zero_cols: ", non_zero_cols
        if tables:
            for name in tables:
                # Get value and save as an array
                vals = tables[name]
                # Generate array of values (FIAT returns [dof, quad_points] transpose to [quad_points, dof])
                value = tabulate_matrix(vals, format)
                code += [(Indent.indent(name), Indent.indent(value))]
                # Tabulate non-zero indices
                if non_zero_cols and self.optimise_level >= 6:
                    code += [Indent.indent(format["comment"]("Array of non-zero columns") )]
                    if name in non_zero_cols:
                        i, cols = non_zero_cols[name]
                        value = format["block"](format["separator"].join(["%d" %c for c in list(cols)]))
                        name_col = format["static const uint declaration"] + format["nonzero columns"](i) + format["array access"](len(cols))
                        code += [(Indent.indent(name_col), value)]  + [""]

                        # Strip the names in the non_zero_cols dictionary for declaration
                        if non_zero_cols and self.optimise_level >= 6:
                            new_name = name.split(" ")[-1].split("[")[0]
                            vals = non_zero_cols[name]
                            del non_zero_cols[name]
                            non_zero_cols[new_name] = vals
                else:
                    code += [""]
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
        return (code, cols_name_map)

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

    def __element_tensor(self, tensor, tensor_number, Indent, format, name_map):
        "Generate loop over primary indices"

        code = []
#        s_set = Set()

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
        elif self.optimise_level == 4:
            values, secondary_loop, t_set, if_statements = values_level_4(indices, vindices, aindices, b0indices,\
                                            bgindices, tensor, tensor_number, weight, format, name_map)
        elif self.optimise_level == 5:
            values, secondary_loop, t_set, if_statements = values_level_5(indices, vindices, aindices, b0indices,\
                                            bgindices, tensor, tensor_number, weight, format, name_map)
        else:
            raise RuntimeError, "Optimisation level not implemented! %d" % self.optimise_level

        value = format_add(values)
        # Count operations to compute entry
        ops_count = reduce_operations.operation_count(value, format)

        # Update ops_count with number of secondary loops
        for sl in secondary_loop:
            ops_count *= sl[2]

        code += [Indent.indent(format["comment"]\
                    ("Number of operations to compute entry =  %d" % (ops_count)))]
        name = ""
        # FIXME: quadrature only support Functionals and Linear and Bilinear forms
        if (irank == 0):

            # Entry is zero because functional is a scalar value
            entry = "0"
            # Generate name
            name =  format["element tensor quad"] + format["array access"](entry)
            code += [Indent.indent(format["comment"]\
                    ("Compute value (tensor/monomial term %d)" % (tensor_number,)))]

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

            # Generate name
            name =  format["element tensor quad"] + format["array access"](entry)
            code += [Indent.indent(format["comment"]\
                    ("Compute block entries (tensor/monomial term %d)" % (tensor_number,)))]

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
        else:
            raise RuntimeError, "Quadrature only support Functionals and Linear and Bilinear forms"

        if self.optimise_level >= 4:
            lines = []
            if not len(values) == len(if_statements):
                raise RuntimeError("The number of if statements and values are not equal.")
            for i in range(len(values)):
                lines.append(if_statements[i])
                lines.append(format["add equal"](name, values[i]))

            if secondary_loop:
                code += generate_loop2(lines, secondary_loop, Indent, format)
            else:
                code += lines

#            if secondary_loop:
#                code += generate_loop(name, value, secondary_loop, Indent, format, format["add equal"])
#            else:
#                code += [format["add equal"](Indent.indent(name), value)]
        else:
            if secondary_loop:
                code += generate_loop(name, value, secondary_loop, Indent, format, format["add equal"])
            else:
                code += [format["add equal"](Indent.indent(name), value)]

        return (code, t_set, ops_count)

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

    def __element_tensor2(self, tensors, group_tensors, Indent, format, psi_name_map, weight_name_map, non_zero_columns, cols_name_map):

        # Prefetch formats to speed up code generation
        format_comment      = format["comment"]
        format_loop         = format["loop"]
        format_ip           = format["integration points"]
        format_block_begin  = format["block begin"]
        format_block_end    = format["block end"]
        format_add          = format["add"]
        format_add_equal    = format["add equal"]
        format_tensor       = format["element tensor quad"]
        format_array_access = format["array access"]


        element_code = []
        trans_set = Set()
#        print "non_zero_columns: ", non_zero_columns
#        print "cols_name_map: ", cols_name_map

        tensor_ops_count = 0
        geo_terms = {}
        for points in group_tensors:
            ip_code = []
            # Loop all quadrature points
            # Create list of tensors for comment
            ts = [group_tensors[points][idims][i] for idims in group_tensors[points]\
                                              for i in range(len(group_tensors[points][idims]))]
            ip_code += [Indent.indent(format["comment"]\
            ("Loop quadrature points (tensor/monomial terms %s)" %(str(tuple(ts)))))]
            ip_code += [Indent.indent(format_loop(format_ip, 0, points))]
            ip_code += [Indent.indent(format_block_begin)]

            # Increase indentation
            Indent.increase()

#            print group_tensors[points]
            # Generate all terms for the given number of quadrature points
            terms, t_set = generate_terms(group_tensors[points], tensors, format, psi_name_map,
                                          weight_name_map, non_zero_columns, cols_name_map, geo_terms, self.optimise_level)
            trans_set = trans_set | t_set
#            print "terms: ", terms
            ip_ops_count = 0
            for prim_loop in terms:
                prim_code = []
                name_dict = terms[prim_loop]

                # Extract inner loops (group those that loop the same indices)
                inner_loops = {}
                for name in name_dict:
#                    print "name: ", name
#                    print name_dict[name]
                    sec_loop_dict = name_dict[name]
                    for sec_loop in sec_loop_dict:
                        # Possibly optimise vals before writing to loop
                        val = format_add(name_dict[name][sec_loop])
                        # Reduce operatoins is optimisation level is 7 or higher 
                        if self.optimise_level >= 7:
                            val = reduce_operations.expand_operations(val, format)
                            if self.optimise_level == 7:
                                val = reduce_operations.reduce_operations(val, format)
                            if self.optimise_level == 8:
                                val = reduce_operations.get_geo_terms(val, geo_terms, self.optimise_level, format)
                            if self.optimise_level == 9:
                                val = reduce_operations.get_geo_terms(val, geo_terms, self.optimise_level, format)
                                val = reduce_operations.reduce_operations(val, format)
                            if self.optimise_level == 10:
                                val = reduce_operations.get_geo_terms(val, geo_terms, self.optimise_level, format)
                                val = reduce_operations.reduce_operations(val, format)
                            if self.optimise_level == 11:
                                val = reduce_operations.get_geo_terms(val, geo_terms, self.optimise_level, format)
                                val = reduce_operations.reduce_operations(val, format)
                        ops_count = reduce_operations.operation_count(val, format)
                        comment = format_comment\
                        ("Number of operations to compute entry = %d" %(ops_count))
                        if not sec_loop in inner_loops:
                            entry =  format_tensor + format_array_access(name)
                            inner_loops[sec_loop] = [comment, format_add_equal(entry, val), ""]
                        else:
                            entry =  format_tensor + format_array_access(name)
                            inner_loops[sec_loop] += [comment, format_add_equal(entry, val), ""]

                        # Update ops count for each of the primary and secondary loop indices
                        for dim in prim_loop[0]:
                            ops_count *= dim
                        for dim in sec_loop[0]:
                            ops_count *= dim
#                        print "sec_loop: ", sec_loop
                        ip_ops_count += ops_count
#                print inner_loops
                lines = []
                for sec_loop in inner_loops:
#                    print "sec loop", sec_loop
                    lines += generate_loop3(inner_loops[sec_loop], sec_loop, Indent, format)
#                keys_vals = get_dict_keys_vals(name_dict, [])
#                for key, val in keys_vals:
#                    name =  format_tensor + format_array_access(key[0])
#                    vals = format_add(val)
#                    line = format_add_equal(name, vals)
#                    lines += generate_loop3([line], key[1], Indent, format)

                prim_code += generate_loop3(lines, prim_loop, Indent, format)

                ip_code += prim_code

            # End the quadrature loop
            # Decrease indentation
            Indent.decrease()
            ip_code += [Indent.indent(format_block_end)]

            comment = format_comment\
                ("Number of operations to compute element tensor for following IP loop = %d" %(ip_ops_count*points) )

            # The operations are done for each IP
            tensor_ops_count += ip_ops_count*points

            element_code += [comment] + ip_code
#        print "geo_terms: ", geo_terms
        # Tabulate goe code
        geo_code = []
        items = geo_terms.items()
#        items = [(v, k) for (k, v) in items]
#        items.sort()
#        items = [(k, v) for (v, k) in items]
        for key, val in items:
            geo_code += [(format["const float declaration"] + geo_terms[key], key)]

        element_code = [format_comment\
        ("Total number of operations to compute element tensor for all IP loops = %d" %(tensor_ops_count) )] + element_code

#            if i + 1 < len(tensors):
#                element_code += [""]

        return (geo_code + element_code, trans_set)





