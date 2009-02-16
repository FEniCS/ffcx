"Utility functions for UFL quadrature code generation"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2009-02-09 -- 2009-02-09"
__copyright__ = "Copyright (C) 2009 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Python modules
from numpy import shape, transpose

# UFL modules
#from ufl.classes import AlgebraOperator, FormArgument
from ufl.classes import *
from ufl.algorithms.analysis import *
from ufl.algorithms.transformations import *

# Utility and optimisation functions for quadraturegenerator
from quadraturegenerator_utils import unique_psi_tables, generate_loop

class BasisTables:
    """This class keeps track of all basis functions, which ones are used and
    it generates the names and tabulates their values."""

    def __init__(self):
        "Initialise class"

        self.element_map = {}
        self.name_map = {}
        self.unique_tables = {}
        self.used_tables = set()
        self.primary_loops = []
        self.functions = {}

    def update(self, tables, optimise_level, format):
        "Update tables"
        self.element_map, self.name_map, self.unique_tables =\
        self.__create_psi_tables(tables, optimise_level, format)
#        self.disp()

    def create_basis_function(self, ufl_basis_function, points, optimise_level, format):
        "Create code for basis functions, and update relevant tables of used basis"

        # FIXME: this needs a lot of work, and it might not even be the best
        # way of doing it
        # Only support test and trial
        indices = {-2: format["first free index"], -1: format["second free index"]}

        element_counter = self.element_map[points][(ufl_basis_function.element(), None)]
        loop_index = indices[ufl_basis_function.count()]
#        print "counter: ", counter
#        print "index: ", index
        name = self.__generate_psi_name(element_counter, None, None, ())
        name, non_zeros = self.name_map[name]
        loop_index_range = shape(self.unique_tables[name])[1]
#        print "name: ", name
        # Append the name to the set of used tables
        self.used_tables.add(name)

        # Change mapping
        mapping = (ufl_basis_function.count(), loop_index, loop_index_range)

        # Create matrix access of basis
        format_ip = format["integration points"]
        if points == 1:
            format_ip = "0"

        return (mapping, name + format["matrix access"](format_ip, loop_index))

    def create_function(self, ufl_function, points, optimise_level, format):
        "Create code for basis functions, and update relevant tables of used basis"

        # FIXME: this needs a lot of work, and it might not even be the best
        # way of doing it
        # Pick first free index of secondary type
        # (could use primary indices, but it's better to avoid confusion)
        loop_index = format["free secondary indices"][0]

        # Get basis name and range
        element_counter = self.element_map[points][(ufl_function.element(), None)]
#        print "counter: ", element_counter
#        print "loop_index: ", loop_index
        basis_name = self.__generate_psi_name(element_counter, None, None, ())
        basis_name, non_zeros = self.name_map[basis_name]
        loop_index_range = shape(self.unique_tables[basis_name])[1]
#        print "basis_name: ", basis_name
        # Add basis name to set of used tables
        self.used_tables.add(basis_name)

        # Add matrix access to basis_name such that we create a unique entry
        # for the expression to compute the function value
        # Create matrix access of basis
        format_ip = format["integration points"]
        if points == 1:
            format_ip = "0"
        basis_name += format["matrix access"](format_ip, loop_index)

        # FIXME: Need to take non-zero mappings, components, restricted and QE elements into account
        coefficient = format["coeff"] + format["matrix access"](str(ufl_function.count()), loop_index)
#        print "basis_name: ", basis_name
#        print "coeff: ", coefficient
        function_expr = format["multiply"]([basis_name, coefficient])

        # Check if the expression to compute the function value is already in
        # the dictionary of used function. If not, generate a new name and add
        function_name = format["function value"] + str(len(self.functions))
        if not function_expr in self.functions:
            self.functions[function_expr] = (function_name, loop_index_range)
        else:
            function_name, index_r = self.functions[function_expr]
            # Check just to make sure
            if not index_r == loop_index_range:
                raise RuntimeError("Index ranges does not match")

        return function_name

    def disp(self):
        print "\nBasisTables, element_map:\n", self.element_map
        print "\nBasisTables, name_map:\n", self.name_map
        print "\nBasisTables, unique_tables:\n", self.unique_tables
        print "\nBasisTables, used_tables:\n", self.used_tables

    def __create_psi_tables(self, tables, optimise_level, format):
        "Create names and maps for tables and non-zero entries if appropriate."

#        print "\nQG-utils, psi_tables:\n", tables

        element_map, flat_tables = self.__flatten_psi_tables(tables)
    #    print "\nQG-utils, psi_tables, flat_tables:\n", flat_tables

        # Outsource call to old function from quadrature_utils
        name_map, unique_tables = unique_psi_tables(flat_tables, optimise_level, format)
    #    name_map, new_tables = unique_psi_tables(flat_tables, 1, format)

    #    print "\nQG-utils, psi_tables, unique_tables:\n", unique_tables
    #    print "\nQG-utils, psi_tables, name_map:\n", name_map

        return (element_map, name_map, unique_tables)

    def __flatten_psi_tables(self, tables):
        "Create a 'flat' dictionary of tables with unique names."

    #    print "\nQG-utils, flatten_psi_tables:\n", tables

        flat_tables = {}
        element_map = {}
        counter = 0
        # Loop quadrature points and get element dictionary {elem: {tables}}
        for point, elem_dict in tables.items():
            element_map[point] = {}
    #        print "\nQG-utils, flatten_tables, points:\n", point
    #        print "\nQG-utils, flatten_tables, elem_dict:\n", elem_dict

            # Loop all elements and get all their tables
            for elem, elem_tables in elem_dict.items():
#                print "\nQG-utils, flatten_tables, elem:\n", elem
    #            print "\nQG-utils, flatten_tables, elem[0].value_rank():\n", elem[0].value_rank()
    #            print "\nQG-utils, flatten_tables, elem_tables:\n", elem_tables
                # If the element value rank != 0, we must loop the components
                # before the derivatives
                # (len(UFLelement.value_shape() == FIATelement.value_rank())
    #            if elem[0].value_rank() != 0:
                element_map[point][elem] = counter
                if len(elem[0].value_shape()) != 0:
                    for num_comp, comp in enumerate(elem_tables):
                        for num_deriv in comp:
                            for derivs, psi_table in num_deriv.items():
    #                            print "\nQG-utils, flatten_tables, derivs:\n", derivs
    #                            print "\nQG-utils, flatten_tables, psi_table:\n", psi_table
                                # Verify shape of basis (can be omitted for speed
                                # if needed I think)
                                if shape(psi_table) != 2 and shape(psi_table)[1] != point:
                                    raise RuntimeError(psi_table, "Something is wrong with this table")

                                name = self.__generate_psi_name(counter, elem[1], num_comp, derivs)
    #                            print "Name: ", name
                                if name in flat_tables:
                                    raise RuntimeError(name, "Name is not unique, something is wrong")
                                flat_tables[name] = transpose(psi_table)
                else:
                    for num_deriv in elem_tables:
                        for derivs, psi_table in num_deriv.items():
    #                        print "\nQG-utils, flatten_tables, derivs:\n", derivs
    #                        print "\nQG-utils, flatten_tables, psi_table:\n", psi_table
                            # Verify shape of basis (can be omitted for speed
                            # if needed I think)
                            if shape(psi_table) != 2 and shape(psi_table)[1] != point:
                                raise RuntimeError(psi_table, "Something is wrong with this table")
                            name = self.__generate_psi_name(counter, elem[1], None, derivs)
    #                        print "Name: ", name
                            if name in flat_tables:
                                raise RuntimeError(name, "Name is not unique, something is wrong")
                            flat_tables[name] = transpose(psi_table)
                counter += 1

        return (element_map, flat_tables)

    def __generate_psi_name(self, counter, restriction, component, derivatives):
        """Generate a name for the psi table of the form:
        FE#_R#_C#_D###, where '#' will be an integer value.

        FE  - is a simple counter to distinguish the various basis, it will be
              assigned in an arbitrary fashion.

        R   - denotes restrictions if applicable, 0 or 1.

        C   - is the component number if any (this does not yet take into account
              tensor valued functions)

        D   - is the number of derivatives in each spatial direction if any. If the
              element is defined in 3D, then D012 means d^3(*)/dydz^2."""

        name = "FE%d" %counter
        if restriction:
            name += "_R%d" % restriction
        if component != None:
            name += "_C%d" % component
        if any(derivatives):
            name += "_D" + "".join([str(d) for d in derivatives])

        return name

def generate_code(integrand, basis_tables, points, optimise_level, Indent, format):
    """Generate code from a UFL integral type.
    This function implements the different optimisation strategies."""

    format_weight       = format["weight"]
    format_scale_factor = format["scale factor"]
    format_add          = format["add"]
    format_add_equal    = format["add equal"]
    format_tensor       = format["element tensor quad"]
    format_array_access = format["array access"]
    format_mult         = format["multiply"]
    format_float        = format["floating point"]
    format_float_decl   = format["float declaration"]
    format_r            = format["free secondary indices"][0]

    # Initialise return values
    code = []
    num_ops = 0
    weights_set = set()
    trans_set = set()

    # We currently only support rank 0, 1 and 2 tensors
    indices = {-2: format["first free index"], -1: format["second free index"]}

    print "\nQG-utils, generate_code, integrand.__repr__():\n", integrand.__repr__()

    loop_code = expand_operations(integrand, basis_tables, points, optimise_level, Indent, format)
    print "loop code: ", loop_code

    # Create code for computing function values, sort after loop ranges first
    functions = basis_tables.functions
    print "FUNC: ", functions
    function_list = {}
    for key, val in functions.items():
        if val[1] in function_list:
            function_list[val[1]].append(key)
        else:
            function_list[val[1]] = [key]
    print "function_list: ", function_list
    # Loop ranges and get list of functions
    for r, func in function_list.items():
        decl = []
        compute = []
        # Loop functions
        for f in func:
            name = functions[f][0]
            decl.append((format_float_decl + name, format_float(0)))
            compute.append(format_add_equal(name, f))
        code += decl + generate_loop(compute, [(format_r, 0, r)], Indent, format)

    # Create weight
    # FIXME: This definitely needs a fix
    weight = format_weight(points)
    if points > 1:
        weight += format["array access"](format["integration points"])

    # Generate entries, multiply by weights and sort after primary loops
    loops = {}
    for key, val in loop_code.items():
        print "Key: ", key

        # Multiply by weight and determinate
        # FIXME: This definitely needs a fix
        value = format["multiply"]([val, weight, format_scale_factor])
        weights_set.add(points)
        trans_set.add(format_scale_factor)

        # FIXME: We only support rank 0, 1 and 2
        entry = ""
        loop = ()
        if len(key) == 0:
            entry = "0"
        elif len(key) == 1:
            key = key[0]
            # Checking if the basis is was a test function
            if key[0] != -2:
                raise RuntimeError("Linear forms must be defined using test functions only")
            # FIXME: Need to consider interior facet integrals
            entry = key[1]
            loop = ((indices[key[0]], 0, key[2]),)
        elif len(key) == 2:
            # Extract test and trial loops in correct order and check if for is legal
            key0, key1 = (0, 0)
            for k in key:
                if not k[0] in indices:
                    raise RuntimeError(k, "Bilinear forms must be defined using test and trial functions (index -2, -1)")
                if k[0] == -2:
                    key0 = k
                else:
                    key1 = k
            # FIXME: Need to consider interior facet integrals
            entry = format_add([format_mult([key0[1], str(key1[2])]), key1[1]])
            loop = ((indices[key0[0]], 0, key0[2]), (indices[key1[0]], 0, key1[2]))
        else:
            raise RuntimeError(key, "Only rank 0, 1 and 2 tensors are currently supported")

        if loop not in loops:
            loops[loop] = [ format_add_equal( format_tensor + format_array_access(entry), value) ]
        else:
            loops[loop].append(format_add_equal( format_tensor + format_array_access(entry), value))

    for loop, lines in loops.items():
        code += generate_loop(lines, loop, Indent, format)

    return (code, num_ops, weights_set, trans_set)

def expand_operations(expr, basis_tables, points, optimise_level, Indent, format):

    print "\nQG-utils, expand_operations, expr.__repr__():\n", expr.__repr__()

    code = {}
    if isinstance(expr, AlgebraOperator):
        print "Algebra: ", expr
        code = apply_algebra(expr, basis_tables, points, optimise_level, Indent, format)
    if isinstance(expr, FormArgument):
        print "FormArgument: ", expr
        code = format_argument(expr, basis_tables, points, optimise_level, Indent, format)

    # FIXME: Make sure that WrapperType only covers containers
    if isinstance(expr, WrapperType):
        print "WrapperType: ", expr
        code = apply_wrapper(expr, basis_tables, points, optimise_level, Indent, format)

    return code

def apply_algebra(expr, basis_tables, points, optimise_level, Indent, format):

    code = {}

    if isinstance(expr, Product):
        format_mult = format["multiply"]
#        print "product: ", expr
#        print "operands: ", expr.operands()
        permute = []
        not_permute = []
        for o in expr.operands():
            print "O: ", o
            expanded = expand_operations(o, basis_tables, points, optimise_level, Indent, format)
            print "expanded: ", expanded

            if len(expanded) > 1 or (expanded and expanded.keys()[0] != ()):
                permute.append(expanded)
            elif expanded:
                not_permute.append(expanded[()])

#        print "permute: ", permute
#        print "not_permute: ", not_permute

        permutations = create_permutations(permute)
#        print "permutations: ", permutations

        if permutations:
            for key in permutations:
                code[key] = format_mult(permutations[key] + not_permute)
        else:
            code[()] = format_mult(not_permute)

    if isinstance(expr, Sum):
        format_add = format["add"]
        format_group = format["grouping"]

        print "sum: ", expr
#        print "operands: ", expr.operands()
#        permute = []
#        not_permute = []
        for o in expr.operands():
            print "O: ", o
            expanded = expand_operations(o, basis_tables, points, optimise_level, Indent, format)
            print "expanded: ", expanded

            if not expanded:
                raise RuntimeError("I didn't expect this")

            # Loop the expanded values and if some entries are the same they
            # can be added, otherwise just dump them in the element tensor
            for key, val in expanded.items():
                if key in code:
                    code[key].append(val)
                else:
                    code[key] = [val]

        # Add sums and group if necessary
        for key, val in code.items():
            if len(val) > 1:
                code[key] = format_group(format_add(val))

    if isinstance(expr, IndexSum):
        format_add = format["add"]
        format_group = format["grouping"]

#        print "index sum: ", expr
#        print "dir(index_sum): ", dir(expr)
        summand, index = expr.operands()
#        print "summand.__repr__(): ", summand.__repr__()
#        print "index.__repr__(): ", index.__repr__()
#        print "index_dims: ", expr.index_dimensions()
#        print "free indices: ", expr.free_indices()
#        print "dir(index): ", dir(index)
#        print "summand.index_dims: ", summand.index_dimensions()
#        print "summand.free indices: ", summand.free_indices()
        # Get index range
        # TODO: What is the #$#%$^%#%^ MultiIndex good for????!!
        # Find better way of getting the index range
        index_range = summand.index_dimensions()[index[0]]
        print "index_range: ", index_range
        new = replace(summand, {index: 0})
        print "NEW: ", new


        # TODO: Does indices always start from 0?
        # This is most likely not the best way of creating code for the
        # IndexSum, but otherwise I'll lose information for the recursive
        # stages such that implicit assumptions must be made.
        # Alternatively, code should be generated for each free index and the
        # sum should just pick the correct indices.
        for i in range(index_range):
            print "I: ", i

        # Just testing to see what this will give me
        test = expand_operations(summand, basis_tables, points, optimise_level, Indent, format)
        print "TEST: ", test

#        permute = []
#        not_permute = []
#        for o in expr.operands():
#            print "O: ", o
#            expanded = expand_operations(o, basis_tables, points, optimise_level, Indent, format)
#            print "expanded: ", expanded

#            if not expanded:
#                raise RuntimeError("I didn't expect this")

#            # Loop the expanded values and if some entries are the same they
#            # can be added, otherwise just dump them in the element tensor
#            for key, val in expanded.items():
#                if key in code:
#                    code[key].append(val)
#                else:
#                    code[key] = [val]

#        # Add sums and group if necessary
#        for key, val in code.items():
#            if len(val) > 1:
#                code[key] = format_group(format_add(val))

    return code

def create_permutations(expr):

#    print "create_permutations, expr: ", expr

    # This is probably not used
    if len(expr) == 0:
        return expr
    # Format keys and values to lists and tuples
    if len(expr) == 1:
        new = {}
        for key, val in expr[0].items():
            if not isinstance(key[0], tuple):
                key = (key,)
            if not isinstance(val, list):
                val = [val]
            new[key] = val

        return new
    # Create permutations of two lists
    # TODO: there could be a cleverer way of changing types of keys and vals
    if len(expr) == 2:
        new = {}
        for key0, val0 in expr[0].items():
            if isinstance(key0[0], tuple):
                key0 = list(key0)
            if not isinstance(key0, list):
                key0 = [key0]
            if not isinstance(val0, list):
                val0 = [val0]
            for key1, val1 in expr[1].items():
                if isinstance(key1[0], tuple):
                    key1 = list(key1)
                if not isinstance(key1, list):
                    key1 = [key1]
                if not isinstance(val1, list):
                    val1 = [val1]
                if tuple(key0 + key1) in new:
                    raise RuntimeError("This is not supposed to happen.")
                new[tuple(key0 + key1)] = val0 + val1

        return new

    # Create permutations by calling this function recursively
    # This is only used for rank > 2 tensors I think
    if len(expr) > 2:
        new = permutations(expr[0:2])
        return permutations(new + expr[2:])

def format_argument(expr, basis_tables, points, optimise_level, Indent, format):

    code = ""
    entry_mapping = ()
    if isinstance(expr, Function):
        print "function: ", expr
#        code = str(expr)
        code = basis_tables.create_function(expr, points, optimise_level, format)

    if isinstance(expr, BasisFunction):
        print "basis function: ", expr
        entry_mapping, code = basis_tables.create_basis_function(expr, points, optimise_level, format)

    return {entry_mapping: code}

def apply_wrapper(expr, basis_tables, points, optimise_level, Indent, format):

    code = {}
    print "\nwrap, dir(expr): ", dir(expr)
    print "wrap, expr.operands(): ", expr.operands()

    # TODO: Is this correct?
    indexed_expr, index = expr.operands()

    print "wrap, indexed_expr: ", indexed_expr
    print "wrap, index: ", index
    print "wrap, expr.free_indices(): ", expr.free_indices()
    print "wrap, expr.index_dimensions(): ", expr.index_dimensions()

    # Get code
    code = expand_operations(indexed_expr, basis_tables, points, optimise_level, Indent, format)

    # Do something to code according to index

    return code








