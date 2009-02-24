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
from ufl.common import *
from ufl.algorithms import *
from ufl.algorithms.ad import *
#from ufl.algorithms.analysis import *
from ufl.algorithms.transformations import *

from ffc.compiler.representation.tensor.multiindex import MultiIndex as FFCMultiIndex
from ffc.compiler.language.tokens import Transform

# Utility and optimisation functions for quadraturegenerator
from quadraturegenerator_utils import unique_psi_tables, generate_loop

class QuadratureTransformer(Transformer):
    "Transform UFL representation to quadrature code"

    def __init__(self, tables, optimise_level, format):

        Transformer.__init__(self)

        # Save format and optimise_level
        self.format = format
        self.optimise_level = optimise_level

        # Create containers and variables
        self.used_tables = set()
        self.primary_loops = []
        self.functions = {}
        self.points = 0

        # Stacks
        self._derivatives = []
        self._index_values = StackDict()
        self._components = StackDict()
        self.trans_set = set()

        self.element_map, self.name_map, self.unique_tables =\
        self.__create_psi_tables(tables, optimise_level, format)

    def reset(self, points):
        # Reset containers
        self.used_tables = set()
        self.primary_loops = []
        self.functions = {}
        self.points = points
        if self._index_values:
            raise RuntimeError("This dictionary is supposed to be empty")
        if self._components:
            raise RuntimeError("This list is supposed to be empty")
        # It should be zero but clear just to be sure
        self._index_values.clear()
        self._components = []

    # Handle the basics just in case, probably not needed?
    def expr(self, o, *operands):
        print "\nVisiting basic Expr:", o.__repr__(), "with operands:"
        print ", ".join(map(str,operands))
        return o

    # Handle the basics just in case, probably not needed?
    def terminal(self, o, *operands):
        print "\nVisiting basic Terminal:", o.__repr__(), "with operands:"
        print ", ".join(map(str,operands))
        return o

    # Algebra
    def sum(self, o, *operands):
        print "\nVisiting Sum:", o.__repr__(), "with operands:"
        print ", ".join(map(str,operands))

        format_group  = self.format["grouping"]
        format_add    = self.format["add"]

        code = {}
        # Loop operands that has to be summend
        for o in operands:
            # If entries does already exist we can add the code, otherwise just
            # dump them in the element tensor
            for key, val in o.items():
                if key in code:
                    code[key].append(val)
                else:
                    code[key] = [val]

        # Add sums and group if necessary
        for key, val in code.items():
            if len(val) > 1:
                code[key] = format_group(format_add(val))
            else:
                code[key] = val[0]

        return code

    def product(self, o, *operands):
        print "\nVisiting Product:", o.__repr__(), "with operands:"
        print ", ".join(map(str,operands))

        format_mult = self.format["multiply"]
        permute = []
        not_permute = []
        # Sort operands in objects that needs permutation and objects that
        # does not
        for op in operands:
            if len(op) > 1 or (op and op.keys()[0] != ()):
                permute.append(op)
            elif op:
                not_permute.append(op[()])
        # Create permutations
        permutations = create_permutations(permute)

#        print "permute: ", permute
#        print "not_permute: ", not_permute
#        print "permutations: ", permutations

        # Create code
        code ={}
        if permutations:
            for key in permutations:
                code[key] = format_mult(permutations[key] + not_permute)
        else:
            code[()] = format_mult(not_permute)

        return code

    def index_sum(self, o):
        print "\nVisiting IndexSum:", o.__repr__(),

        summand, index = o.operands()
        print "\nindex.__repr__(): ", index.__repr__()
        print "\nindex[0].__repr__(): ", index[0].__repr__()
        # Just a few safety checks
        if not isinstance(index, MultiIndex) and len(index) == 1:
            raise RuntimeError(index, "Expecting 1 MultiIndex")

        format_group = self.format["grouping"]
        format_add = self.format["add"]
        tmp = 0
        code = {}
        for i in range(o.dimension()):
            self._index_values.push(index[0], i)
            tmp = self.visit(summand)
#            print "index_sum, tmp: ", tmp
            self._index_values.pop()

            # FIXME: remove this?
            if not tmp:
                raise RuntimeError("I didn't expect this")

            # If entries does already exist we can add the code, otherwise just
            # dump them in the element tensor
            for key, val in tmp.items():
                if key in code:
                    code[key].append(val)
                else:
                    code[key] = [val]

        # Add sums and group if necessary
        for key, val in code.items():
            if len(val) > 1:
                code[key] = format_group(format_add(val))
            else:
                code[key] = val

        return code

    def indexed(self, o):
        print "\nVisiting Indexed:", o.__repr__(),

        indexed_expr, index = o.operands()
        print "wrap, indexed_expr: ", indexed_expr
        print "wrap, index.__repr__(): ", index.__repr__()

        # Loop multi indices and create components
        for indx in index:
            self._components.append(indx)
            # If index is not present in index_values, create it.
            # (It means that it is a Fixed index??!)
            if not indx in self._index_values:
                if not isinstance(indx, FixedIndex):
                    raise RuntimeError(indx, "Index must be Fixed for Indexed to add it to index_values")
                self._index_values.push(indx, indx._value)

        # Visit expression subtrees and generate code
        code = self.visit(indexed_expr)

        # Loop multi indices and delete components
        for indx in index:
            self._components.pop()
            if isinstance(indx, FixedIndex):

                self._index_values.pop()

        return code

    def component_tensor(self, o):
        print "\nVisiting ComponentTensor:", o.__repr__(),

        indexed_expr, index = o.operands()
        print "wrap, indexed_expr: ", indexed_expr
        print "wrap, index.__repr__(): ", index.__repr__()

        if not len(self._components) == len(index):
            raise RuntimeError("The number of known components must be equal to the number of components of the ComponentTensor for this to work.")

        # Save copy of components to let parent delete them again
        old_components = self._components[:]
        self._components = []

        # Loop multi indices and map index values
        for i, indx in enumerate(index):
            self._index_values.push(indx, self._index_values[old_components[i]])

        # Visit expression subtrees and generate code
        code = self.visit(indexed_expr)

        # Loop multi indices and delete index values
        for indx in index:
            self._index_values.pop()

        # Set components equal to old components
        self._components = old_components[:]

        return code

    def spatial_derivative(self, o):
        print "\nVisiting SpatialDerivative:", o.__repr__(),

        indexed_expr, index = o.operands()
        print "wrap, indexed_expr: ", indexed_expr
        print "wrap, index: ", index.__repr__()
        if not isinstance(index, MultiIndex) and len(index) == 1:
            raise RuntimeError(index, "Expecting 1 MultiIndex")

        # Get the direction that we need the derivative for
        direction = None
        if isinstance(index[0], FixedIndex):
            direction = index[0]._value
        else:
            direction = self._index_values[index[0]]

        # Append the derivative
        self._derivatives.append(direction)

        # Visit children
        code = self.visit(indexed_expr)
        print "spatial_derivative, code: ", code

        return code

    # FormArguments
#    def basis_function(self, o, *operands, component=[]):
    def basis_function(self, o, *operands):
        print "\nVisiting BasisFunction:", o.__repr__()

        # Just checking that we don't get any operands
        if operands:
            raise RuntimeError(operands, "Didn't expect any operands for BasisFunction")

#        print "self._components: ", self._components
        if len(self._components) > 1:
            raise RuntimeError(self._components, "Currently only supports 1 component value (tensor valued basis not supported)")

#        print "self._derivatives: ", self._derivatives
        # Create aux. info
        # FIXME: restriction not handled yet
        restriction = None
        component = None
        derivatives = ()
        # Handle derivatives and components
        if self._derivatives:
            derivatives = self._derivatives[:]
        if self._components:
            component = self._index_values[self._components[0]]

        print "\nDerivatives: ", derivatives

        # Create mapping and code for basis function
        basis = self.create_basis_function(o, restriction, component, derivatives)

        # Reset spatial derivatives
        # FIXME: (should this be handled by SpatialDerivative)
        self._derivatives = []

        return basis

    def function(self, o, *operands):
        print "\nVisiting Function:", o.__repr__()

        # Just checking that we don't get any operands
        if operands:
            raise RuntimeError(operands, "Didn't expect any operands for BasisFunction")

#        print "self._components: ", self._components
        if len(self._components) > 1:
            raise RuntimeError(self._components, "Currently only supports 1 component value (tensor valued basis not supported)")

#        print "self._derivatives: ", self._derivatives
        # Create aux. info
        # FIXME: restriction not handled yet
        restriction = None
        component = None
        derivatives = ()
        # Handle derivatives and components
        if self._derivatives:
            derivatives = self._derivatives[:]
        if self._components:
            component = self._index_values[self._components[0]]

        print "\nDerivatives: ", derivatives

        # Create code for basis function
        code = self.create_function(o, restriction, component, derivatives)

        # Reset spatial derivatives
        # FIXME: (should this be handled by SpatialDerivative)
        self._derivatives = []

        return {(): code}


    def create_basis_function(self, ufl_basis_function, restriction, component, derivatives):
        "Create code for basis functions, and update relevant tables of used basis"

        format_ip            = self.format["integration points"]
        format_matrix_access = self.format["matrix access"]
        format_group         = self.format["grouping"]
        format_add           = self.format["add"]
        format_mult          = self.format["multiply"]
        format_transform     = self.format["transform"]

        # Only support test and trial functions
        # TODO: Verify that test and trial functions will ALWAYS be rearranged to 0 and 1
        indices = {-2: self.format["first free index"],
                   -1: self.format["second free index"],
                    0: self.format["first free index"],
                    1: self.format["second free index"]}

        element_counter = self.element_map[self.points][(ufl_basis_function.element(), restriction)]
        if not ufl_basis_function.count() in indices:
            raise RuntimeError(ufl_basis_function, "Currently, BasisFunction index must be either -2, -1, 0 or 1")

        loop_index = indices[ufl_basis_function.count()]
#        print "counter: ", counter
#        print "index: ", index

        code = {}
        # Generate FFC multi index for derivatives
#        print "dims: ", [range(ufl_basis_function.element().cell().d)]*sum(derivatives)
        geo_dim = ufl_basis_function.element().cell().d
        multiindices = FFCMultiIndex([range(geo_dim)]*len(derivatives)).indices
#        print "multiindices: ", multiindices
        # Loop derivatives and get multi indices
        for multi in multiindices:
            deriv = [multi.count(i) for i in range(geo_dim)]
            print "multi: ", multi
            if not any(deriv):
                deriv = []
            print "deriv: ", deriv

            name = self.__generate_psi_name(element_counter, restriction, component, deriv)
#        print "name_map: ", self.name_map
            name, non_zeros = self.name_map[name]
            loop_index_range = shape(self.unique_tables[name])[1]
##        print "name: ", name
            # Append the name to the set of used tables
            self.used_tables.add(name)

            # Change mapping
            mapping = ((ufl_basis_function.count(), loop_index, loop_index_range),)

            # Create matrix access of basis
            if self.points == 1:
                format_ip = "0"
            basis = name + format_matrix_access(format_ip, loop_index)

            # Add transformation if supported and needed
            transforms = []
            for i, direction in enumerate(derivatives):
                ref = multi[i]
                if ufl_basis_function.element().family() != "Lagrange":
                    raise RuntimeError(ufl_basis_function.element().family(), "Only derivatives of Lagrange elements is currently supported")
                t = format_transform(Transform.JINV, ref, direction, restriction)
                self.trans_set.add(t)
                transforms.append(t)

            if mapping in code:
                code[mapping].append(format_mult(transforms + [basis]))
            else:
                code[mapping] = [format_mult(transforms + [basis])]

        # Add sums and group if necessary
        for key, val in code.items():
            if len(val) > 1:
                code[key] = format_group(format_add(val))
            else:
                code[key] = val[0]

        return code

    def create_function(self, ufl_function, restriction, component, derivatives):
        "Create code for basis functions, and update relevant tables of used basis"

        format_ip            = self.format["integration points"]
        format_matrix_access = self.format["matrix access"]
        format_group         = self.format["grouping"]
        format_add           = self.format["add"]
        format_mult          = self.format["multiply"]
        format_transform     = self.format["transform"]
        format_coeff         = self.format["coeff"]
        format_F             = self.format["function value"]

        # FIXME: this needs a lot of work, and it might not even be the best
        # way of doing it
        # Pick first free index of secondary type
        # (could use primary indices, but it's better to avoid confusion)
        loop_index = self.format["free secondary indices"][0]

        # Get basis name and range
        element_counter = self.element_map[self.points][(ufl_function.element(), restriction)]
#        print "counter: ", element_counter
#        print "loop_index: ", loop_index

        code = []

        # Generate FFC multi index for derivatives
#        print "dims: ", [range(ufl_basis_function.element().cell().d)]*sum(derivatives)
        geo_dim = ufl_function.element().cell().d
        multiindices = FFCMultiIndex([range(geo_dim)]*len(derivatives)).indices

        for multi in multiindices:
            deriv = [multi.count(i) for i in range(geo_dim)]
            print "multi: ", multi
            if not any(deriv):
                deriv = []
            print "deriv: ", deriv

            basis_name = self.__generate_psi_name(element_counter, restriction, component, deriv)
            basis_name, non_zeros = self.name_map[basis_name]
            loop_index_range = shape(self.unique_tables[basis_name])[1]
            print "basis_name: ", basis_name
            # Add basis name to set of used tables
            self.used_tables.add(basis_name)

            # Add matrix access to basis_name such that we create a unique entry
            # for the expression to compute the function value
            # Create matrix access of basis
            if self.points == 1:
                format_ip = "0"
            basis_name += format_matrix_access(format_ip, loop_index)

            # FIXME: Need to take non-zero mappings, components, restricted and QE elements into account
            coefficient = format_coeff + format_matrix_access(str(ufl_function.count()), loop_index)

            function_expr = format_mult([basis_name, coefficient])

            # Add transformation if supported and needed
            transforms = []
            for i, direction in enumerate(derivatives):
                ref = multi[i]
                if ufl_function.element().family() != "Lagrange":
                    raise RuntimeError(ufl_function.element().family(), "Only derivatives of Lagrange elements is currently supported")
                t = format_transform(Transform.JINV, ref, direction, restriction)
                self.trans_set.add(t)
                transforms.append(t)
            function_expr = format_mult(transforms + [function_expr])

            # Check if the expression to compute the function value is already in
            # the dictionary of used function. If not, generate a new name and add
            function_name = format_F + str(len(self.functions))
            if not function_expr in self.functions:
                self.functions[function_expr] = (function_name, loop_index_range)
            else:
                function_name, index_r = self.functions[function_expr]
                # Check just to make sure
                if not index_r == loop_index_range:
                    raise RuntimeError("Index ranges does not match")
            code.append(function_name)

        if len(code) > 1:
            code = format_group(format_add(code))
        else:
            code = code[0]

        return code

    def disp(self):
        print "\nQuadratureTransformer, element_map:\n", self.element_map
        print "\nQuadratureTransformer, name_map:\n", self.name_map
        print "\nQuadratureTransformer, unique_tables:\n", self.unique_tables
        print "\nQuadratureTransformer, used_tables:\n", self.used_tables

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

def generate_code(integrand, transformer, points, optimise_level, Indent, format):
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

    print "\nQG, Using Transformer"
    # Expand all derivatives
    print "Integrand: ", integrand
    print "Integrand: ", expand_derivatives(integrand)

    transformer.reset(points)
    loop_code = transformer.visit(integrand)
    trans_set = transformer.trans_set
    print "loop_code: ", loop_code

    # TODO: Verify that test and trial functions will ALWAYS be rearranged to 0 and 1
    indices = {-2: format["first free index"], -1: format["second free index"],
                0: format["first free index"],  1: format["second free index"]}

    print "\nQG-utils, generate_code, integrand.__repr__():\n", integrand.__repr__()

    # Create code for computing function values, sort after loop ranges first
    functions = transformer.functions
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

        # Multiply by weight and determinant
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
            # Checking if the basis was a test function
            # TODO: Make sure test function indices are always rearranged to 0
            if key[0] != -2 and key[0] != 0:
                raise RuntimeError("Linear forms must be defined using test functions only")
            # FIXME: Need to consider interior facet integrals
            entry = key[1]
            loop = ((indices[key[0]], 0, key[2]),)
        elif len(key) == 2:
            # Extract test and trial loops in correct order and check if for is legal
            key0, key1 = (0, 0)
            for k in key:
                if not k[0] in indices:
                    raise RuntimeError(k, "Bilinear forms must be defined using test and trial functions (index -2, -1, 0, 1)")
                if k[0] == -2 or k[0] == 0:
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



