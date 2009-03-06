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
#from ufl.algorithms import *
#from ufl.algorithms.ad import expand_derivatives
from ufl.algorithms.analysis import extract_sub_elements
from ufl.algorithms.transformations import *
from ufl.algorithms.printing import tree_format

from ffc.compiler.representation.tensor.multiindex import MultiIndex as FFCMultiIndex
from ffc.compiler.language.tokens import Transform
from ffc.compiler.language.restriction import *

# Utility and optimisation functions for quadraturegenerator
from quadraturegenerator_utils import generate_loop, unique_tables, get_ones, contains_zeros
from reduce_operations import operation_count, expand_operations, reduce_operations

class QuadratureTransformer(Transformer):
    "Transform UFL representation to quadrature code"

    def __init__(self, form_representation, domain_type, optimise_options, format):

        Transformer.__init__(self)

        # Save format, optimise_options, weights and fiat_elements_map
        self.format = format
        self.optimise_options = optimise_options
        self.quadrature_weights = form_representation.quadrature_weights[domain_type]
        self.fiat_elements_map = form_representation.fiat_elements_map

        # Create containers and variables
        self.used_psi_tables = set()
        self.used_weights = set()
        self.trans_set = set()
        self.functions = {}
        self.function_count = 0
        self.geo_dim = 0
        self.points = 0
        self.facet0 = None
        self.facet1 = None
        self.restriction = None

        # Stacks
        self._derivatives = []
        self._index_values = StackDict()
        self._components = StackDict()
        self.trans_set = set()
        self.element_map, self.name_map, self.unique_tables =\
        self.__create_psi_tables(form_representation.psi_tables[domain_type])

    def update_facets(self, facet0, facet1):
        self.facet0 = facet0
        self.facet1 = facet1
        # Reset functions and count everytime we generate a new case of facets
        self.functions = {}
        self.function_count = 0

    def update_points(self, points):
        self.points = points
        # Reset functions everytime we move to a new quadrature loop
        # But not the functions count.
        self.functions = {}

    def reset(self):
        # Reset containers
        self.used_psi_tables = set()
        self.used_weights = set()
        self.trans_set = set()
        self.functions = {}
        self.function_count = 0
        self.geo_dim = 0
        self.points = 0
        self.facet0 = None
        self.facet1 = None
        if self._index_values:
            raise RuntimeError("This dictionary is supposed to be empty")
        if self._components:
            raise RuntimeError("This list is supposed to be empty")
        # It should be zero but clear just to be sure
        self._index_values.clear()
        self._components = []

    def disp(self):
        print "\n\n **** Displaying QuadratureTransformer ****"
        print "\nQuadratureTransformer, element_map:\n", self.element_map
        print "\nQuadratureTransformer, name_map:\n", self.name_map
        print "\nQuadratureTransformer, unique_tables:\n", self.unique_tables
        print "\nQuadratureTransformer, used_psi_tables:\n", self.used_psi_tables
        print "\nQuadratureTransformer, used_weights:\n", self.used_weights

    # Handle the basics just in case, probably not needed?
    def expr(self, o, *operands):
        print "\n\nVisiting basic Expr:", o.__repr__(), "with operands:"
        print ", ".join(map(str,operands))
#        return {}
        # FIXME: unsafe switch back on
        return o

    # Handle the basics just in case, probably not needed?
    def terminal(self, o, *operands):
        print "\n\nVisiting basic Terminal:", o.__repr__(), "with operands:"
        print ", ".join(map(str,operands))
#        return {}
        # FIXME: unsafe switch back on
        return o

    # -------------------------------------------------------------------------
    # Things which are not supported (geometry.py, expr)
    # -------------------------------------------------------------------------
    def facet_normal(self, o):
        print "\n\nVisiting FacetNormal:", o.__repr__()
        raise RuntimeError("FacetNormal is not supported (yet), use a VectorElement(Discontinuous Lagrange, 'cell_type', 0) instead.")

    def wrapper_type(self, o):
        print "\n\nVisiting WrapperType:", o.__repr__()
        raise RuntimeError(o.__repr__(), "This WrapperType is not supported (yet).")

    # -------------------------------------------------------------------------
    # Constant values (constantvalue.py)
    # -------------------------------------------------------------------------
    def float_value(self, o, *operands):
#        print "\n\nVisiting FloatValue:", o.__repr__()

        # FIXME: Might be needed because it can be IndexAnnotated?
        if operands:
            raise RuntimeError((o, operands), "Did not expect any operands for FloatValue")

        # TODO: Handle value < 0

        return {():self.format["floating point"](o.value())}

    def int_value(self, o, *operands):
#        print "\n\nVisiting IntValue:", o.__repr__()

        # FIXME: Might be needed because it can be IndexAnnotated?
        if operands:
            raise RuntimeError((o, operands), "Did not expect any operands for IntValue")

        # TODO: Handle value < 0

        return {():self.format["floating point"](o.value())}

    def identity(self, o, *operands):
        print "\n\nVisiting Identity:", o.__repr__()
        raise RuntimeError("Identity should have been expanded!!")

    # -------------------------------------------------------------------------
    # Algebra (algebra.py)
    # -------------------------------------------------------------------------
    def sum(self, o, *operands):
        debug("Visiting Sum: " + o.__repr__() + "\noperands: " + "\n".join(map(str, operands)))

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
#        print "\n\nVisiting Product:", o.__repr__(), "with operands:"
#        print ", ".join(map(str,operands))

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
            for key, val in permutations.items():
                # Sort key in order to create a unique key
                l = list(key)
                l.sort()
                code[tuple(l)] = format_mult(val + not_permute)
        else:
            code[()] = format_mult(not_permute)

        return code

    def division(self, o, *operands):
#        print "\n\nVisiting Division:", o.__repr__(), "with operands:"
#        print ", ".join(map(str,operands))

        format_div = self.format["division"]

        if len(operands) != 2:
            raise RuntimeError(operands, "Expected exactly two operands (numerator and denominator) ")

        numerator_code, denominator_code = operands
#        print "numerator: ", numerator_code
#        print "denominator: ", denominator_code

        # TODO: Are these safety checks needed?
        if not () in denominator_code and len(denominator_code) != 1:
            raise RuntimeError(denominator_code, "Only support function type denominator")

        denominator = denominator_code.pop(())

        for key, val in numerator_code.items():
            numerator_code[key] = val + format_div + denominator

        return numerator_code


    def power(self, o):
#        print "\n\nVisiting Power:", o.__repr__()

        # Get base and exponent
        base, expo = o.operands()
#        print "expo: ", expo

        # Visit base
        base_code = self.visit(base)
#        print "base: ", base_code

        # TODO: Are these safety checks needed?
        if not () in base_code and len(base_code) != 1:
            raise RuntimeError(base_code, "Only support function type base")

        val = base_code.pop(())

        # Multiply val by self expo number of times
        return {(): self.format["power"](val, expo.value())}

    def abs(self, o, *operands):
        print "\n\nVisiting Abs:", o.__repr__(), "with operands:"
        print ", ".join(map(str,operands))

        # TODO: Are these safety checks needed?
        if len(operands) != 1 and not () in operands[0] and len(operands[0]) != 1:
            raise RuntimeError(operands, "Abs expects one operand of function type")

        # Take absolute value of operand
        operand = operands[0]
        for key, val in operand.items():
            operand[key] = self.format["absolute value"](val)

        return operand

    # -------------------------------------------------------------------------
    # MathFunctions (mathfunctions.py)
    # -------------------------------------------------------------------------
    def sqrt(self, o, *operands):
        print "\n\nVisiting Sqrt:", o.__repr__(), "with operands:"
        print ", ".join(map(str,operands))

        # TODO: Are these safety checks needed?
        if len(operands) != 1 and not () in operands[0] and len(operands[0]) != 1:
            raise RuntimeError(operands, "Sqrt expects one operand of function type")

        # Take absolute value of operand
        operand = operands[0]
        for key, val in operand.items():
            operand[key] = self.format["sqrt"](val)

        return operand

    def exp(self, o, *operands):
        print "\n\nVisiting Exp:", o.__repr__(), "with operands:"
        print ", ".join(map(str,operands))

        # TODO: Are these safety checks needed?
        if len(operands) != 1 and not () in operands[0] and len(operands[0]) != 1:
            raise RuntimeError(operands, "Exp expects one operand of function type")

        # Take absolute value of operand
        operand = operands[0]
        for key, val in operand.items():
            operand[key] = self.format["exp"](val)

        return operand

    def ln(self, o, *operands):
        print "\n\nVisiting Ln:", o.__repr__(), "with operands:"
        print ", ".join(map(str,operands))

        # TODO: Are these safety checks needed?
        if len(operands) != 1 and not () in operands[0] and len(operands[0]) != 1:
            raise RuntimeError(operands, "Ln expects one operand of function type")

        # Take absolute value of operand
        operand = operands[0]
        for key, val in operand.items():
            operand[key] = self.format["ln"](val)

        return operand

    def cos(self, o, *operands):
        print "\n\nVisiting Cos:", o.__repr__(), "with operands:"
        print ", ".join(map(str,operands))

        # TODO: Are these safety checks needed?
        if len(operands) != 1 and not () in operands[0] and len(operands[0]) != 1:
            raise RuntimeError(operands, "Cos expects one operand of function type")

        # Take absolute value of operand
        operand = operands[0]
        for key, val in operand.items():
            operand[key] = self.format["cos"](val)

        return operand

    def sin(self, o, *operands):
        print "\n\nVisiting Sin:", o.__repr__(), "with operands:"
        print ", ".join(map(str,operands))

        # TODO: Are these safety checks needed?
        if len(operands) != 1 and not () in operands[0] and len(operands[0]) != 1:
            raise RuntimeError(operands, "Sin expects one operand of function type")

        # Take absolute value of operand
        operand = operands[0]
        for key, val in operand.items():
            operand[key] = self.format["sin"](val)

        return operand

    # -------------------------------------------------------------------------
    # Retriction (restriction.py)
    # -------------------------------------------------------------------------
    def positive_restricted(self, o):
#        print "\n\nVisiting PositiveRestricted:", o.__repr__(),

        restricted_expr = o.operands()
#        print "\noperands", restricted_expr
        if len(restricted_expr) != 1:
            raise RuntimeError(restricted_expr, "Only expected one operand for restriction")
 
        # Just visit the first operand, there should only be one
        # FIXME: Need to handle restriction
        self.restriction = Restriction.PLUS
        code = self.visit(restricted_expr[0])
        # Reset restriction
        # TODO: Is this really necessary?
        self.restriction = None
        return code

    def negative_restricted(self, o):
#        print "\n\nVisiting NegativeRestricted:", o.__repr__(),

        # Just get the first operand, there should only be one
        restricted_expr = o.operands()
#        print "\noperands", restricted_expr
        if len(restricted_expr) != 1:
            raise RuntimeError(restricted_expr, "Only expected one operand for restriction")
 
        # Just visit the first operand, there should only be one
        # FIXME: Need to handle restriction
        self.restriction = Restriction.MINUS
        code = self.visit(restricted_expr[0])
        # Reset restriction
        # TODO: Is this really necessary?
        self.restriction = None
        return code

    # -------------------------------------------------------------------------
    # From indexed.py, indexsum.py and tensors.py
    # -------------------------------------------------------------------------
    def index_sum(self, o):
        debug("\n\nVisiting IndexSum:" + o.__repr__())

        raise RuntimeError("IndexSum should have been expanded")
        summand, index = o.operands()
#        print "\nindex.__repr__(): ", index.__repr__()
#        print "\nindex[0].__repr__(): ", index[0].__repr__()
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
        debug("\n\nVisiting Indexed:" + o.__repr__())

        indexed_expr, index = o.operands()

#        print "\nwrap, indexed_expr: ", indexed_expr.__repr__()
#        print "\nwrap, index.__repr__(): ", index.__repr__()

        # Save copy of components to let parent delete them again
        old_components = self._components[:]
        self._components = []

        # Loop multi indices and create components
        for indx in index:
            self._components.append(indx)
            # If index is not present in index_values, create it.
            # (It means that it is a Fixed index??!)
            if not indx in self._index_values:
                if not isinstance(indx, FixedIndex):
                    raise RuntimeError(indx, "Index must be Fixed for Indexed to add it to index_values")
                self._index_values.push(indx, indx._value)

#        print "\ncomponents: ", self._components
#        print "\nindex_values: ", self._index_values

        # Visit expression subtrees and generate code
        code = self.visit(indexed_expr)

        # Loop multi indices and delete components
        for indx in index:
#            self._components.pop()
            if isinstance(indx, FixedIndex):
                self._index_values.pop()

        # Set components equal to old components
        self._components = old_components[:]

        return code

    def component_tensor(self, o):
        debug("\n\nVisiting ComponentTensor:" + o.__repr__())

        indexed_expr, index = o.operands()
#        print "wrap, indexed_expr: ", indexed_expr
#        print "wrap, index.__repr__(): ", index.__repr__()

        if not len(self._components) == len(index):
            raise RuntimeError("The number of known components must be equal to the number of components of the ComponentTensor for this to work.")

        # Save copy of components to let parent delete them again
        old_components = self._components[:]
        self._components = []

        # Loop multi indices and map index values
        for i, indx in enumerate(index):
            self._index_values.push(indx, self._index_values[old_components[i]])

#        print "\nCTcomponents: ", self._components
#        print "\nCTindex_values: ", self._index_values

        # Visit expression subtrees and generate code
        code = self.visit(indexed_expr)

        # Loop multi indices and delete index values
        for indx in index:
            self._index_values.pop()

        # Set components equal to old components
        self._components = old_components[:]

        return code

    def list_tensor(self, o):
        print "\n\nVisiting ListTensor:", o.__repr__(),

        raise RuntimeError("ListTensors should have been expanded!")
        print "\nLTcomponents: ", self._components
        print "\nLTindex_values: ", self._index_values

        # A list tensor only has one component
        if len(self._components) != 1:
            raise RuntimeError(self._components, "ListTensor can only have one component")

        expr = o.operands()[self._index_values[self._components[0]]]

        return self.visit(expr)

    # -------------------------------------------------------------------------
    # Derivatives (differentiation.py)
    # -------------------------------------------------------------------------
    def spatial_derivative(self, o):
#        print "\n\nVisiting SpatialDerivative:", o.__repr__(),

        indexed_expr, index = o.operands()
#        print "wrap, indexed_expr: ", indexed_expr
#        print "wrap, index: ", index.__repr__()
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
#        print "spatial_derivative, code: ", code

        return code

    # -------------------------------------------------------------------------
    # BasisFunction, Function and Constants (basisfunction.py, function.py)
    # -------------------------------------------------------------------------
    def basis_function(self, o, *operands):
#        print "\n\nVisiting BasisFunction:", o.__repr__()

        # Just checking that we don't get any operands
        if operands:
            raise RuntimeError(operands, "Didn't expect any operands for BasisFunction")

#        print "self._components: ", self._components
        if len(self._components) > 1:
            raise RuntimeError(self._components, "Currently only supports 1 component value (tensor valued basis not supported)")

#        print "self._derivatives: ", self._derivatives
        # Create aux. info
        component = None
        derivatives = ()
        # Handle derivatives and components
        if self._derivatives:
            derivatives = self._derivatives[:]
        if self._components:
            component = self._index_values[self._components[0]]

#        print "\nDerivatives: ", derivatives

        # Create mapping and code for basis function
        basis = self.create_basis_function(o, component, derivatives)

        # Reset spatial derivatives
        # FIXME: (should this be handled by SpatialDerivative)
        self._derivatives = []

        return basis

    def constant(self, o, *operands):
#        print "\n\nVisiting Constant:", o.__repr__()
        # Just checking that we don't get any operands
        if operands:
            raise RuntimeError(operands, "Didn't expect any operands for Constant")

        if len(self._components) > 0:
            raise RuntimeError(self._components, "Constant does not expect components")

        # FIXME: restriction not handled yet
        restriction = None
        coefficient = self.format["coeff"] + self.format["matrix access"](str(o.count()), 0)
        return {():coefficient}

    def vector_constant(self, o, *operands):
#        print "\n\nVisiting VectorConstant:", o.__repr__()
        # Just checking that we don't get any operands
        if operands:
            raise RuntimeError(operands, "Didn't expect any operands for VectorConstant")

        if len(self._components) != 1:
            raise RuntimeError(self._components, "VectorConstant only expext 1 component")

        # FIXME: restriction not handled yet
        restriction = None

        # We get one component
        component = self._index_values[self._components[0]]
        coefficient = self.format["coeff"] + self.format["matrix access"](str(o.count()), component)
        return {():coefficient}

    def function(self, o, *operands):
#        print "\n\nVisiting Function:", o.__repr__()

        # Just checking that we don't get any operands
        if operands:
            raise RuntimeError(operands, "Didn't expect any operands for Function")

#        print "self._components: ", self._components
        if len(self._components) > 1:
            raise RuntimeError(self._components, "Currently only supports 1 component value (tensor valued functions not supported)")

#        print "self._derivatives: ", self._derivatives
        # Create aux. info
        component = None
        derivatives = ()
        # Handle derivatives and components
        if self._derivatives:
            derivatives = self._derivatives[:]
        if self._components:
            component = self._index_values[self._components[0]]

#        print "\nDerivatives: ", derivatives

        # Create code for basis function
        code = self.create_function(o, component, derivatives)

        # Reset spatial derivatives
        # FIXME: (should this be handled by SpatialDerivative)
        self._derivatives = []

        return {(): code}

    # -------------------------------------------------------------------------
    # Helper functions for BasisFunction and Function)
    # -------------------------------------------------------------------------
    def create_basis_function(self, ufl_basis_function, component, derivatives):
        "Create code for basis functions, and update relevant tables of used basis"

        format_ip            = self.format["integration points"]
        format_matrix_access = self.format["matrix access"]
        format_array_access  = self.format["array access"]
        format_group         = self.format["grouping"]
        format_add           = self.format["add"]
        format_mult          = self.format["multiply"]
        format_transform     = self.format["transform"]
        format_nzc           = self.format["nonzero columns"]

        # Only support test and trial functions
        # TODO: Verify that test and trial functions will ALWAYS be rearranged to 0 and 1
        indices = {-2: self.format["first free index"],
                   -1: self.format["second free index"],
                    0: self.format["first free index"],
                    1: self.format["second free index"]}

        # Check that we have a basis function
        if not ufl_basis_function.count() in indices:
            raise RuntimeError(ufl_basis_function, "Currently, BasisFunction index must be either -2, -1, 0 or 1")

        # Check that we don't take derivatives of QuadratureElements
        # FIXME: We just raise an exception now, but should we just return 0?
        # UFL will apply Dx(f_e*f_qe, 0) which should result in f_e.dx(0)*f_qe*dx
        # instead of an error?
        if derivatives and any(e.family() == "Quadrature" for e in extract_sub_elements(ufl_basis_function.element())):
            raise RuntimeError(ufl_basis_function, "Derivatives of Quadrature elements are not supported")

        # Handle restriction through facet
        facet = {Restriction.PLUS: self.facet0, Restriction.MINUS: self.facet1, None: self.facet0}[self.restriction]

        # Get element counter and loop index
        element_counter = self.element_map[self.points][ufl_basis_function.element()]
        loop_index = indices[ufl_basis_function.count()]

        # Create basis access, we never need to map the entry in the basis
        # table since we will either loop the entire space dimension or the
        # non-zeros
        if self.points == 1:
            format_ip = "0"
        basis_access = format_matrix_access(format_ip, loop_index)


        # Get the FIAT element and offset by element space dimension in case of
        # negative restriction
        fiat_element = self.fiat_elements_map[ufl_basis_function.element()]
        space_dim = fiat_element.space_dimension()
        offset = {Restriction.PLUS: "", Restriction.MINUS: str(space_dim), None: ""}[self.restriction]
        # If we have a restricted function multiply space_dim by two
        if self.restriction == Restriction.PLUS or self.restriction == Restriction.MINUS:
            space_dim *= 2

#        print "counter: ", counter
#        print "index: ", index

        code = {}
        # Generate FFC multi index for derivatives

        # Set geo_dim
        # TODO: All terms REALLY have to be defined on cell with the same
        # geometrical dimension so only do this once and exclude the check?
        geo_dim = ufl_basis_function.element().cell().d
        if self.geo_dim:
            if geo_dim != self.geo_dim:
                raise RuntimeError(geo_dim, "All terms must be defined on cells with the same geometrical dimension")
        else:
            self.geo_dim = geo_dim

        # Generate FFC multi index for derivatives
        multiindices = FFCMultiIndex([range(geo_dim)]*len(derivatives)).indices
#        print "multiindices: ", multiindices

        # Loop derivatives and get multi indices
        for multi in multiindices:
            deriv = [multi.count(i) for i in range(geo_dim)]
#            print "multi: ", multi
            if not any(deriv):
                deriv = []
#            print "deriv: ", deriv

            # TODO: Handle zeros and ones info
            name = self.__generate_psi_name(element_counter, facet, component, deriv)
            name, non_zeros, zeros, ones = self.name_map[name]
#            print "\nname: ", name
#            print "\nnon_zeros: ", non_zeros
#            print "\nzeros: ", zeros
#            print "\nones: ", ones
            loop_index_range = shape(self.unique_tables[name])[1]

            # Append the name to the set of used tables and create
            # matrix access
            # TODO: Handle this more elegantly such that all terms involving this
            # zero factor is removed
            basis = "0"
            if not (zeros and self.optimise_options["ignore zero tables"] == 1):
                self.used_psi_tables.add(name)
                basis = name + basis_access

            # Create the correct mapping of the basis function into the local
            # element tensor
            basis_map = loop_index
            if non_zeros:
                basis_map = format_nzc(non_zeros[0]) + format_array_access(basis_map)
            if offset:
                basis_map = format_group(format_add([basis_map, offset]))

            # Create mapping (index, map, loop_range, space_dim)
            # Example dx and ds: (0, j, 3, 3)
            # Example dS: (0, (j + 3), 3, 6), 6=2*space_dim
            # Example dS optimised: (0, (nz2[j] + 3), 2, 6), 6=2*space_dim
            mapping = ((ufl_basis_function.count(), basis_map, loop_index_range, space_dim),)

            # Add transformation if supported and needed
            transforms = []
            for i, direction in enumerate(derivatives):
                ref = multi[i]
                # FIXME: Handle other element types too
                if ufl_basis_function.element().family() not in ["Lagrange", "Discontinuous Lagrange"]:
                    if ufl_basis_function.element().family() == "Mixed":
                        # Check that current sub components only contain supported elements
#                        print "Component: ", component
#                        print "basis: ", basis
                        if not all(e.family() in ["Lagrange", "Discontinuous Lagrange", "Mixed"] for e in extract_sub_elements(ufl_basis_function.element())):
                            raise RuntimeError(ufl_basis_function.element().family(), "Only derivatives of Lagrange elements is currently supported")
                    else:
                        raise RuntimeError(ufl_basis_function.element().family(), "Only derivatives of Lagrange elements is currently supported")
                t = format_transform(Transform.JINV, ref, direction, self.restriction)
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

    def create_function(self, ufl_function, component, derivatives):
        "Create code for basis functions, and update relevant tables of used basis"

        format_ip            = self.format["integration points"]
        format_matrix_access = self.format["matrix access"]
        format_array_access = self.format["array access"]
        format_group         = self.format["grouping"]
        format_add           = self.format["add"]
        format_mult          = self.format["multiply"]
        format_transform     = self.format["transform"]
        format_coeff         = self.format["coeff"]
        format_F             = self.format["function value"]
        format_nzc           = self.format["nonzero columns"]

        # Check that we don't take derivatives of QuadratureElements
        # FIXME: We just raise an exception now, but should we just return 0?
        # UFL will apply Dx(f_e*f_qe, 0) which should result in f_e.dx(0)*f_qe*dx
        # instead of an error?
        # Is this even correct? What if we just want the first component of a
        # mixed element which is not a quadrature element?
        quad_element = any(e.family() == "Quadrature" for e in extract_sub_elements(ufl_function.element()))
        if derivatives and quad_element:
            raise RuntimeError(ufl_function, "Derivatives of Quadrature elements are not supported")

        # Pick first free index of secondary type
        # (could use primary indices, but it's better to avoid confusion)
        loop_index = self.format["free secondary indices"][0]

        # Create basis access, we never need to map the entry in the basis
        # table since we will either loop the entire space dimension or the
        # non-zeros
        if self.points == 1:
            format_ip = "0"
        basis_access = format_matrix_access(format_ip, loop_index)

        # Handle restriction through facet
        facet = {Restriction.PLUS: self.facet0, Restriction.MINUS: self.facet1, None: self.facet0}[self.restriction]

        # Get the element counter
        element_counter = self.element_map[self.points][ufl_function.element()]

        # Get the FIAT element and offset by element space dimension in case of
        # negative restriction
        fiat_element = self.fiat_elements_map[ufl_function.element()]
        offset = {Restriction.PLUS: "", Restriction.MINUS: str(fiat_element.space_dimension()), None: ""}[self.restriction]

        code = []

#        print "dims: ", [range(ufl_basis_function.element().cell().d)]*sum(derivatives)
        # Set geo_dim
        # TODO: All terms REALLY have to be defined on cell with the same
        # geometrical dimension so only do this once and exclude the check?
        geo_dim = ufl_function.element().cell().d
        if self.geo_dim:
            if geo_dim != self.geo_dim:
                raise RuntimeError(geo_dim, "All terms must be defined on cells with the same geometrical dimension")
        else:
            self.geo_dim = geo_dim

        # Generate FFC multi index for derivatives
        multiindices = FFCMultiIndex([range(geo_dim)]*len(derivatives)).indices
        for multi in multiindices:
            deriv = [multi.count(i) for i in range(geo_dim)]
#            print "multi: ", multi
            if not any(deriv):
                deriv = []
#            print "deriv: ", deriv

            # TODO: Handle zeros info
            basis_name = self.__generate_psi_name(element_counter, facet, component, deriv)
            basis_name, non_zeros, zeros, ones = self.name_map[basis_name]
#            print "\nbasis_name: ", basis_name
#            print "\nnon_zeros: ", non_zeros
#            print "\nzeros: ", zeros
#            print "\nones: ", ones
            # If all basis are zero we just return "0"
            # TODO: Handle this more elegantly such that all terms involving this
            # zero factor is removed
            if zeros and self.optimise_options["ignore zero tables"]:
                continue

            # Get the index range of the loop index
            loop_index_range = shape(self.unique_tables[basis_name])[1]
#            print "\nloop index range: ", loop_index_range

            # Set default coefficient access
            coefficient_access = loop_index

            # If the loop index range is one we can look up the first component
            # in the coefficient array. If we only have ones we don't need the basis
            if self.optimise_options["ignore ones"] > 0 and loop_index_range == 1 and ones:
                coefficient_access = "0"
                basis_name = ""
            else:
                # Add basis name to set of used tables and add
                # matrix access
                self.used_psi_tables.add(basis_name)
                basis_name += basis_access

            # If we have a quadrature element we can use the ip number to look
            # up the value directly. Need to add offset in case of components
            if quad_element:
                quad_offset = 0
                if component:
                    for i in range(component):
                        quad_offset += fiat_element.sub_element(i).space_dimension()
                if quad_offset:
                    coefficient_access = format_add([format_ip, str(quad_offset)])
                else:
                    coefficient_access = format_ip
            # If we have non zero column mapping but only one value just pick it
            if non_zeros and coefficient_access == "0":
                coefficient_access = str(non_zeros[1][0])
            elif non_zeros:
                coefficient_access = format_nzc(non_zeros[0]) + format_array_access(coefficient_access)
            if offset:
                coefficient_access = format_add([coefficient_access, offset])

            # Try to evaluate coefficient access ("3 + 2" --> "5")
            try:
                coefficient_access = str(eval(coefficient_access))
            except:
                pass

            coefficient = format_coeff + format_matrix_access(str(ufl_function.count()), coefficient_access)
            function_expr = coefficient
            if basis_name:
                function_expr = format_mult([basis_name, coefficient])

            # Add transformation if supported and needed
            transforms = []
            for i, direction in enumerate(derivatives):
                ref = multi[i]
                # FIXME: Handle other element types too
                if ufl_function.element().family() not in ["Lagrange", "Discontinuous Lagrange"]:
                    if ufl_function.element().family() == "Mixed":
                        # Check that current sub components only contain supported elements
#                        print "Component: ", component
#                        print "basis: ", basis
                        if not all(e.family() in ["Lagrange", "Discontinuous Lagrange", "Mixed"] for e in extract_sub_elements(ufl_function.element())):
                            raise RuntimeError(ufl_function.element().family(), "Only derivatives of Lagrange elements is currently supported")
                    else:
                        raise RuntimeError(ufl_function.element().family(), "Only derivatives of Lagrange elements is currently supported")
                # Create transform and add to set of used transformations
                t = format_transform(Transform.JINV, ref, direction, self.restriction)
                self.trans_set.add(t)
                transforms.append(t)

            # If we have a quadrature element (or if basis was deleted) we
            # don't need the basis
            if quad_element or not basis_name:
                function_name = coefficient
            else:
                # Check if the expression to compute the function value is already in
                # the dictionary of used function. If not, generate a new name and add
                function_name = format_F + str(self.function_count)
                if not function_expr in self.functions:
                    self.functions[function_expr] = (function_name, loop_index_range)
                    # Increase count
                    self.function_count += 1
                else:
                    function_name, index_r = self.functions[function_expr]
                    # Check just to make sure
                    if not index_r == loop_index_range:
                        raise RuntimeError("Index ranges does not match")

            # Multiply function value by the transformations and add to code
            code.append(format_mult(transforms + [function_name]))

        if not code:
            return "0"
        elif len(code) > 1:
            code = format_group(format_add(code))
        else:
            code = code[0]

        return code

    def __create_psi_tables(self, tables):
        "Create names and maps for tables and non-zero entries if appropriate."

#        print "\nQG-utils, psi_tables:\n", tables

        element_map, flat_tables = self.__flatten_psi_tables(tables)
    #    print "\nQG-utils, psi_tables, flat_tables:\n", flat_tables

        # Outsource call to old function from quadrature_utils
        name_map, unique_tables = self.__unique_psi_tables(flat_tables)
    #    name_map, new_tables = unique_psi_tables(flat_tables, 1, format)

#        print "\nQG-utils, psi_tables, unique_tables:\n", unique_tables
#        print "\nQG-utils, psi_tables, name_map:\n", name_map

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
#            print "\nQG-utils, flatten_tables, points:\n", point
#            print "\nQG-utils, flatten_tables, elem_dict:\n", elem_dict

            # Loop all elements and get all their tables
            for elem, facet_tables in elem_dict.items():
#                print "\nQG-utils, flatten_tables, elem:\n", elem
#                print "\nQG-utils, flatten_tables, facet_tables:\n", facet_tables
                # If the element value rank != 0, we must loop the components
                # before the derivatives
                # (len(UFLelement.value_shape() == FIATelement.value_rank())
                element_map[point][elem] = counter
                for facet, elem_tables in facet_tables.items():
                    if len(elem.value_shape()) != 0:
                        for num_comp, comp in enumerate(elem_tables):
                            for num_deriv in comp:
                                for derivs, psi_table in num_deriv.items():
#                                    print "\nQG-utils, flatten_tables, derivs:\n", derivs
#                                    print "\nQG-utils, flatten_tables, psi_table:\n", psi_table
                                    # Verify shape of basis (can be omitted for speed
                                    # if needed I think)
                                    if shape(psi_table) != 2 and shape(psi_table)[1] != point:
                                        raise RuntimeError(psi_table, "Something is wrong with this table")

                                    name = self.__generate_psi_name(counter, facet, num_comp, derivs)
#                                    print "Name: ", name
                                    if name in flat_tables:
                                        raise RuntimeError(name, "Name is not unique, something is wrong")
                                    flat_tables[name] = transpose(psi_table)
                    else:
                        for num_deriv in elem_tables:
                            for derivs, psi_table in num_deriv.items():
#                                print "\nQG-utils, flatten_tables, derivs:\n", derivs
#                                print "\nQG-utils, flatten_tables, psi_table:\n", psi_table
                                # Verify shape of basis (can be omitted for speed
                                # if needed I think)
                                if shape(psi_table) != 2 and shape(psi_table)[1] != point:
                                    raise RuntimeError(psi_table, "Something is wrong with this table")
                                name = self.__generate_psi_name(counter, facet, None, derivs)
#                                print "Name: ", name
                                if name in flat_tables:
                                    raise RuntimeError(name, "Name is not unique, something is wrong")
                                flat_tables[name] = transpose(psi_table)
                counter += 1
#        raise RuntimeError
        return (element_map, flat_tables)

    def __generate_psi_name(self, counter, facet, component, derivatives):
        """Generate a name for the psi table of the form:
        FE#_f#_C#_D###, where '#' will be an integer value.

        FE  - is a simple counter to distinguish the various basis, it will be
              assigned in an arbitrary fashion.

        f   - denotes facets if applicable, range(element.num_facets()).

        C   - is the component number if any (this does not yet take into account
              tensor valued functions)

        D   - is the number of derivatives in each spatial direction if any. If the
              element is defined in 3D, then D012 means d^3(*)/dydz^2."""

        name = "FE%d" %counter
        if facet != None:
            name += "_f%d" % facet
        if component != None:
            name += "_C%d" % component
        if any(derivatives):
            name += "_D" + "".join([str(d) for d in derivatives])

        return name

    def __unique_psi_tables(self, tables):
        "Determine if some tensors have the same tables (and same names)"

        # Get formats
        format_epsilon = self.format["epsilon"]

        # Get unique tables
        name_map, inverse_name_map = unique_tables(tables, format_epsilon)

#        print "\nname_map: ", name_map
#        print "\ninv_name_map: ", inverse_name_map
#        print "\ntables: ", tables

        # Get names of tables with all ones
        names_ones = get_ones(tables, format_epsilon)

        # Set values to zero if they are lower than threshold
        for name in tables:
            # Get values
            vals = tables[name]
            for r in range(shape(vals)[0]):
                for c in range(shape(vals)[1]):
                    if abs(vals[r][c]) < format_epsilon:
                        vals[r][c] = 0
            tables[name] = vals

        # Extract the column numbers that are non-zero
        # (only for optimisations higher than 0)
        i = 0
        non_zero_columns = {}
        if self.optimise_options["non zero columns"]:
            for name in tables:
                # Get values
                vals = tables[name]

                # Use the first row as reference
                non_zeros = list(vals[0].nonzero()[0])

                # If all columns in the first row are non zero, there's no point
                # in continuing
                if len(non_zeros) == shape(vals)[1]:
                    continue

                # If we only have one row (IP) we just need the nonzero columns
                if shape(vals)[0] == 1:
                    if list(non_zeros):
                        non_zeros.sort()
                        non_zero_columns[name] = (i, non_zeros)

                        # Compress values
                        tables[name] = vals[:, non_zeros]
                        i += 1

                # Check if the remaining rows are nonzero in the same positions, else expand
                else:
                    for j in range(shape(vals)[0] - 1):
                        # All rows must have the same non-zero columns
                        # for the optimization to work (at this stage)
                        new_non_zeros = list(vals[j+1].nonzero()[0])
                        if non_zeros != new_non_zeros:
                            non_zeros = non_zeros + [c for c in new_non_zeros if not c in non_zeros]
                            # If this results in all columns being non-zero, continue.
                            if len(non_zeros) == shape(vals)[1]:
                                continue

                    # Only add nonzeros if it implies a reduction of columns
                    if len(non_zeros) != shape(vals)[1]:
                        if list(non_zeros):
                            non_zeros.sort()
                            non_zero_columns[name] = (i, non_zeros)

                            # Compress values
                            tables[name] = vals[:, non_zeros]
                            i += 1

        # Check if we have some zeros in the tables
        names_zeros = self.__contains_zeros(tables, format_epsilon)

        # Add non-zero column info to inverse_name_map
        # (so we only need to pass around one name_map to code generating functions)
        for name in inverse_name_map:
            if inverse_name_map[name] in non_zero_columns:
                nzc = non_zero_columns[inverse_name_map[name]]
                zero = inverse_name_map[name] in names_zeros
                ones = inverse_name_map[name] in names_ones
#                print "zero: ", zero
#                print "ones: ", ones
                inverse_name_map[name] = [inverse_name_map[name], nzc, zero, ones]
            else:
                zero = inverse_name_map[name] in names_zeros
                ones = inverse_name_map[name] in names_ones
#                print "zero: ", zero
#                print "ones: ", ones
                inverse_name_map[name] = [inverse_name_map[name], (), zero, ones]

        # If we found non zero columns we might be able to reduce number of tables
        # further if optimise level is higher than 0
        if non_zero_columns:

            # Try reducing the tables. This is possible if some tables have become
            # identical as a consequence of compressing the tables
            nm, inv_nm = unique_tables(tables, format_epsilon)

            # Update name maps
            for name in inverse_name_map:
                if inverse_name_map[name][0] in inv_nm:
                    inverse_name_map[name][0] = inv_nm[inverse_name_map[name][0]]
            for name in nm:
                maps = nm[name]
                for m in maps:
                    if not name in name_map:
                        name_map[name] = []
                    if m in name_map:
                        name_map[name] += name_map[m] + [m]
                        del name_map[m]
                    else:
                        name_map[name].append(m)

            # Exclude tables with all ones
            names = get_ones(tables, format_epsilon)
            # Because these tables now contain ones as a consequence of compression
            # we still need to consider the non-zero columns when looking up values
            # in coefficient arrays. The psi entries can however we neglected and we
            # don't need to tabulate the values
            for name in names:
                if name in name_map:
                    maps = name_map[name]
                    for m in maps:
                        inverse_name_map[m][3] = True
                if name in inverse_name_map:
                        inverse_name_map[name][3] = True

        # Write protect info
        for name in inverse_name_map:
            inverse_name_map[name] = tuple(inverse_name_map[name])

        return (inverse_name_map, tables)

    def __contains_zeros(self, tables, format_epsilon):
        "Checks if any tables contains all zeros"

        names = []
        for name in tables:
            vals = tables[name]
            zero = True
            for r in range(shape(vals)[0]):
                if not zero:
                    break
                for c in range(shape(vals)[1]):
                    # If just one value is different from zero, break loops
                    if abs(vals[r][c]) > format_epsilon:
                        zero = False
                        break

            if zero:
                print "\n*** Warning: this table only contains zeros. This is not critical,"
                print "but it might slow down the runtime performance of your code!"
                print "Do you take derivatives of a constant?"
                names.append(name)
        return names

def generate_code(integrand, transformer, Indent, format):
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
    format_comment      = format["comment"]
    format_F            = format["function value"]

    # Initialise return values
    code = []
    num_ops = 0

#    print "\nQG, Using Transformer"
    # Apply basic expansions
    # TODO: Figure out if there is a 'correct' order of doing this
    # In form.form_data().form, which we should be using, coefficients have
    # been mapped and derivatives expande. So it should be enough to just
    # expand_indices and purge_list_tensors
#    print "Integrand: ", integrand
    new_integrand = expand_indices(integrand)
#    print "Integrand: ", new_integrand
    new_integrand = purge_list_tensors(new_integrand)
#    print "Integrand: ", new_integrand
#    print "Expanded integrand\n", tree_format(new_integrand)

    loop_code = transformer.visit(new_integrand)
#    print "loop_code: ", loop_code

    # TODO: Verify that test and trial functions will ALWAYS be rearranged to 0 and 1
    indices = {-2: format["first free index"], -1: format["second free index"],
                0: format["first free index"],  1: format["second free index"]}

#    print "\nQG-utils, generate_code, integrand.__repr__():\n", integrand.__repr__()

    # Create code for computing function values, sort after loop ranges first
    functions = transformer.functions
    function_list = {}
    for key, val in functions.items():
        if val[1] in function_list:
            function_list[val[1]].append(key)
        else:
            function_list[val[1]] = [key]
#    print "function_list: ", function_list

    # Create the function declarations, we know that the code generator numbers
    # functions from 0 to n.
    if transformer.function_count:
        code += ["", format_comment("Function declarations")]
    for function_number in range(transformer.function_count):
        code.append((format_float_decl + format_F + str(function_number), format_float(0)))

    # Loop ranges and get list of functions
    for loop_range, list_of_functions in function_list.items():
        function_expr = {}
        function_numbers = []
        func_ops = 0
        # Loop functions
        for function in list_of_functions:
            # Get name and number
            name = functions[function][0]
            number = int(name.strip(format_F))
            # TODO: This check can be removed for speed later
            if number in function_numbers:
                raise RuntimeError("This is definitely not supposed to happen!")
            function_numbers.append(number)
            # Get number of operations to compute entry and add to function
            # operations count
            f_ops = operation_count(function, format) + 1
            func_ops += f_ops
            function_expr[number] = format_add_equal(name, function)

        # Multiply number of operations by the range of the loop index and add
        # number of operations to compute function values to total count
        func_ops *= loop_range
        func_ops_comment = ["", format_comment("Total number of operations to compute function values = %d" % func_ops)]
        num_ops += func_ops

        # Sort the functions according to name and create loop to compute the
        # function values
        function_numbers.sort()
        lines = []
        for number in function_numbers:
            lines.append(function_expr[number])
        code += func_ops_comment + generate_loop(lines, [(format_r, 0, loop_range)], Indent, format)

    # Create weight
    # FIXME: This definitely needs a fix
    weight = format_weight(transformer.points)
    if transformer.points > 1:
        weight += format["array access"](format["integration points"])

    # Generate entries, multiply by weights and sort after primary loops
    loops = {}
    for key, val in loop_code.items():
#        print "Key: ", key

        # If value was zero continue
        if val == None:
            continue
        # Multiply by weight and determinant
        # FIXME: This definitely needs a fix
        value = format_mult([val, weight, format_scale_factor])
        transformer.used_weights.add(transformer.points)
        transformer.trans_set.add(format_scale_factor)

        # Use old operation reduction
        if transformer.optimise_options["simplify expressions"] == 2:
            value = expand_operations(value, format)
            value = reduce_operations(value, format)

        # Compute number of operations to compute entry and create comment
        # (add 1 because of += in assignment)
        entry_ops = operation_count(value, format) + 1
        entry_ops_comment = format_comment("Number of operations to compute entry = %d" % entry_ops)
        prim_ops = entry_ops

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
                raise RuntimeError(key, "Linear forms must be defined using test functions only")

            index_j, entry, range_j, space_dim_j = key
            loop = ((indices[index_j], 0, range_j),)
            # Multiply number of operations to compute entries by range of loop
            prim_ops *= range_j
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
            index_j, entry_j, range_j, space_dim_j = key0
            index_k, entry_k, range_k, space_dim_k = key1

            entry = format_add([format_mult([entry_j, str(space_dim_k)]), entry_k])
            loop = ((indices[index_j], 0, range_j), (indices[index_k], 0, range_k))

            # Multiply number of operations to compute entries by range of loops
            prim_ops *= range_j*range_k
        else:
            raise RuntimeError(key, "Only rank 0, 1 and 2 tensors are currently supported")

        # Generate the code line for the entry
        entry_code = format_add_equal( format_tensor + format_array_access(entry), value)

        if loop not in loops:
            loops[loop] = [prim_ops, [entry_ops_comment, entry_code]]
        else:
            loops[loop][0] += prim_ops
            loops[loop][1] += [entry_ops_comment, entry_code]

    for loop, ops_lines in loops.items():
        ops, lines = ops_lines
        # Add number of operations for current loop to total count
        num_ops += ops
        code += ["", format_comment("Number of operations for primary indices = %d" % ops)]
        code += generate_loop(lines, loop, Indent, format)

    return (code, num_ops)

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



