"QuadratureTransformerBase, a common class for quadrature transformers to translate UFL expressions."

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@gmail.com)"
__date__ = "2009-10-13"
__copyright__ = "Copyright (C) 2009-2010 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-02-08

# Python modules.
from itertools import izip
import time
from numpy import shape

# UFL Classes.
from ufl.classes import MultiIndex, FixedIndex, Index
from ufl.common import StackDict, Stack

# UFL Algorithms.
from ufl.algorithms import propagate_restrictions, Transformer, tree_format, strip_variables

# FFC modules.
from ffc.log import ffc_assert, error, info
from ffc.fiatinterface import create_element
from ffc.mixedelement import MixedElement
from ffc.cpp import format

# FFC tensor modules.
from ffc.tensor.multiindex import MultiIndex as FFCMultiIndex

# Utility and optimisation functions for quadraturegenerator.
from quadratureutils import create_psi_tables
from symbolics import generate_aux_constants
from symbolics import BASIS, IP, GEO, CONST

class QuadratureTransformerBase(Transformer):
#class QuadratureTransformerBase(ReuseTransformer):
    "Transform UFL representation to quadrature code."

    def __init__(self, psi_tables, optimise_parameters):

        Transformer.__init__(self)

        # Save optimise_parameters, weights and fiat_elements_map.
        self.optimise_parameters = optimise_parameters

        # Create containers and variables.
        self.used_psi_tables = set()
        self.psi_tables_map = {}
        self.used_weights = set()
        self.used_nzcs = set()
        self.geo_consts = {}
        self.ip_consts = {}
        self.trans_set = set()
        self.functions = {}
        self.function_count = 0
        self.geo_dim = 0
        self.points = 0
        self.facet0 = None
        self.facet1 = None
        self.restriction = None

        # Stacks.
        self._derivatives = []
        self._index2value = StackDict()
        self._components = Stack()
        self.element_map, self.name_map, self.unique_tables =\
              create_psi_tables(psi_tables, self.optimise_parameters)

        # Cache.
        self.argument_cache = {}
        self.function_cache = {}

    def update_facets(self, facet0, facet1):
        self.facet0 = facet0
        self.facet1 = facet1
#        # Reset functions and count everytime we generate a new case of facets.
#        self.functions = {}
#        self.function_count = 0

#        # Reset cache
#        self.argument_cache = {}
#        self.function_cache = {}

    def update_points(self, points):
        self.points = points
        # Reset functions everytime we move to a new quadrature loop
        self.functions = {}
        self.function_count = 0

        # Reset cache
        self.argument_cache = {}
        self.function_cache = {}

    def disp(self):
        print "\n\n **** Displaying QuadratureTransformer ****"
        print "\nQuadratureTransformer, element_map:\n", self.element_map
        print "\nQuadratureTransformer, name_map:\n", self.name_map
        print "\nQuadratureTransformer, unique_tables:\n", self.unique_tables
        print "\nQuadratureTransformer, used_psi_tables:\n", self.used_psi_tables
        print "\nQuadratureTransformer, psi_tables_map:\n", self.psi_tables_map
        print "\nQuadratureTransformer, used_weights:\n", self.used_weights
        print "\nQuadratureTransformer, geo_consts:\n", self.geo_consts

    def component(self):
        "Return current component tuple."
        if len(self._components):
            return self._components.peek()
        return ()

    def derivatives(self):
        "Return all derivatives tuple."
        if len(self._derivatives):
            return tuple(self._derivatives[:])
        return ()

    # -------------------------------------------------------------------------
    # Start handling UFL classes.
    # -------------------------------------------------------------------------
    # Nothing in expr.py is handled. Can only handle children of these clases.
    def expr(self, o):
        print "\n\nVisiting basic Expr:", repr(o), "with operands:"
        error("This expression is not handled: ", repr(o))

    # Nothing in terminal.py is handled. Can only handle children of these clases.
    def terminal(self, o):
        print "\n\nVisiting basic Terminal:", repr(o), "with operands:"
        error("This terminal is not handled: ", repr(o))

    # -------------------------------------------------------------------------
    # Things which should not be here (after expansion etc.) from:
    # algebra.py, differentiation.py, finiteelement.py,
    # form.py, geometry.py, indexing.py, integral.py, tensoralgebra.py, variable.py.
    # -------------------------------------------------------------------------
    def algebra_operator(self, o, *operands):
        print "\n\nVisiting AlgebraOperator: ", repr(o)
        error("This type of AlgebraOperator should have been expanded!!" + repr(o))

    def derivative(self, o, *operands):
        print "\n\nVisiting Derivative: ", repr(o)
        error("All derivatives apart from SpatialDerivative should have been expanded!!")

    def finite_element_base(self, o, *operands):
        print "\n\nVisiting FiniteElementBase: ", repr(o)
        error("FiniteElements must be member of a Argument or Coefficient!!")

    def form(self, o, *operands):
        print "\n\nVisiting Form: ", repr(o)
        error("The transformer only work on a Form integrand, not the Form itself!!")

    def space(self, o):
        print "\n\nVisiting Space: ", repr(o)
        error("A Space should not be present in the integrand.")

    def cell(self, o):
        print "\n\nVisiting Cell: ", repr(o)
        error("A Cell should not be present in the integrand.")

    def index_base(self, o):
        print "\n\nVisiting IndexBase: ", repr(o)
        error("Indices should not be floating around freely in the integrand!!")

    def integral(self, o):
        print "\n\nVisiting Integral: ", repr(o)
        error("Integral should not be present in the integrand!!")

    def measure(self, o):
        print "\n\nVisiting Measure: ", repr(o)
        error("Measure should not be present in the integrand!!")

    def compound_tensor_operator(self, o):
        print "\n\nVisiting CompoundTensorOperator: ", repr(o)
        error("CompoundTensorOperator should have been expanded.")

    def label(self, o):
        print "\n\nVisiting Label: ", repr(o)
        error("What is a Lable doing in the integrand?")

    # -------------------------------------------------------------------------
    # Things which are not supported yet, from:
    # condition.py, constantvalue.py, function.py, geometry.py, lifting.py,
    # mathfunctions.py, restriction.py
    # -------------------------------------------------------------------------
    def condition(self, o):
        print "\n\nVisiting Condition:", repr(o)
        error("Condition is not supported (yet).")

    def conditional(self, o):
        print "\n\nVisiting Condition:", repr(o)
        error("Conditional is not supported (yet).")

    def constant_value(self, o):
        print "\n\nVisiting ConstantValue:", repr(o)
        error("This type of ConstantValue is not supported (yet).")

    def index_annotated(self, o):
        print "\n\nVisiting IndexAnnotated:", repr(o)
        error("Only child classes of IndexAnnotated is supported.")

    def constant_base(self, o):
        print "\n\nVisiting ConstantBase:", repr(o)
        error("This type of ConstantBase is not supported (yet).")

    def geometric_quantity(self, o):
        print "\n\nVisiting GeometricQuantity:", repr(o)
        error("This type of GeometricQuantity is not supported (yet).")

    def spatial_coordinate(self, o):
        print "\n\nVisiting SpatialCoordinate:", repr(o)
        error("SpatialCoordinate is not supported (yet).")

    def lifting_result(self, o):
        print "\n\nVisiting LiftingResult:", repr(o)
        error("LiftingResult (and children) is not supported (yet).")

    def terminal_operator(self, o):
        print "\n\nVisiting TerminalOperator:", repr(o)
        error("TerminalOperator (LiftingOperator and LiftingFunction) is not supported (yet).")

    def math_function(self, o):
        print "\n\nVisiting MathFunction:", repr(o)
        error("This MathFunction is not supported (yet).")

    def restricted(self, o):
        print "\n\nVisiting Restricted:", repr(o)
        error("This type of Restricted is not supported (only positive and negative are currently supported).")

    # -------------------------------------------------------------------------
    # Handlers that should be implemented by child classes.
    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------
    # AlgebraOperators (algebra.py).
    # -------------------------------------------------------------------------
    def sum(self, o, *operands):
        print "\n\nVisiting Sum: ", repr(o)
        error("This object should be implemented by the child class.")

    def product(self, o, *operands):
        print "\n\nVisiting Product: ", repr(o)
        error("This object should be implemented by the child class.")

    def division(self, o, *operands):
        print "\n\nVisiting Division: ", repr(o)
        error("This object should be implemented by the child class.")

    def power(self, o):
        print "\n\nVisiting Power: ", repr(o)
        error("This object should be implemented by the child class.")

    def abs(self, o, *operands):
        print "\n\nVisiting Abs: ", repr(o)
        error("This object should be implemented by the child class.")

    # -------------------------------------------------------------------------
    # FacetNormal (geometry.py).
    # -------------------------------------------------------------------------
    def facet_normal(self, o,  *operands):
        print "\n\nVisiting FacetNormal: ", repr(o)
        error("This object should be implemented by the child class.")

    # -------------------------------------------------------------------------
    # Things that can be handled by the base class.
    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------
    # Argument (basisfunction.py).
    # -------------------------------------------------------------------------
    def argument(self, o, *operands):
        #print("\nVisiting Argument:" + repr(o))

        # Just checking that we don't get any operands.
        ffc_assert(not operands, "Didn't expect any operands for Argument: " + repr(operands))

        # Create aux. info.
        components = self.component()
        derivatives = self.derivatives()

        # Check if basis is already in cache
        basis = self.argument_cache.get((o, components, derivatives, self.restriction), None)
        # FIXME: Why does using a code dict from cache make the expression manipulations blow (MemoryError) up later?
        if basis is not None and not self.optimise_parameters["simplify expressions"]:
#        if basis is not None:
            return basis

        # Get auxiliary variables to generate basis
        component, local_comp, local_offset, ffc_element, quad_element, \
        transformation, multiindices = self._get_auxiliary_variables(o, components, derivatives)

        # Create mapping and code for basis function and add to dict.
        basis = self.create_argument(o, derivatives, component, local_comp,
                  local_offset, ffc_element, transformation, multiindices)

        self.argument_cache[(o, components, derivatives, self.restriction)] = basis

        return basis

    # -------------------------------------------------------------------------
    # Constant values (constantvalue.py).
    # -------------------------------------------------------------------------
    def identity(self, o):
        #print "\n\nVisiting Identity: ", repr(o)

        # Get components
        components = self.component()

        # Safety checks.
        ffc_assert(not o.operands(), "Didn't expect any operands for Identity: " + repr(o.operands()))
        ffc_assert(len(components) == 2, "Identity expect exactly two component indices: " + repr(components))

        # Only return a value if i==j
        if components[0] == components[1]:
            return self._format_scalar_value(1.0)
        return self._format_scalar_value(None)

    def scalar_value(self, o, *operands):
        "ScalarValue covers IntValue and FloatValue"
        #print "\n\nVisiting ScalarValue: ", repr(o)

        # FIXME: Might be needed because it can be IndexAnnotated?
        ffc_assert(not operands, "Did not expect any operands for ScalarValue: " + repr((o, operands)))
        return self._format_scalar_value(o.value())

    def zero(self, o, *operands):
        #print "\n\nVisiting Zero:", repr(o)
        # FIXME: Might be needed because it can be IndexAnnotated?
        ffc_assert(not operands, "Did not expect any operands for Zero: " + repr((o, operands)))
        return self._format_scalar_value(None)

    # -------------------------------------------------------------------------
    # SpatialDerivative (differentiation.py).
    # -------------------------------------------------------------------------
    def spatial_derivative(self, o):
        #print("\n\nVisiting SpatialDerivative: " + repr(o))

        # Get expression and index
        derivative_expr, index = o.operands()

        # Get direction of derivative and check that we only get one return index
        der = self.visit(index)
        ffc_assert(len(der) == 1, "SpatialDerivative: expected only one direction index. " + repr(der))

        # Add direction to list of derivatives
        self._derivatives.append(der[0])

        # Visit children to generate the derivative code.
        code = self.visit(derivative_expr)

        # Remove the direction from list of derivatives
        self._derivatives.pop()
        return code

    # -------------------------------------------------------------------------
    # Coefficient and Constants (function.py).
    # -------------------------------------------------------------------------
    def coefficient(self, o, *operands):
        #print("\nVisiting Coefficient: " + repr(o))

        # Safety check.
        ffc_assert(not operands, "Didn't expect any operands for Coefficient: " + repr(operands))

        # Create aux. info.
        components = self.component()
        derivatives = self.derivatives()

        # Check if function is already in cache
        function_code = self.function_cache.get((o, components, derivatives, self.restriction), None)
        # FIXME: Why does using a code dict from cache make the expression manipulations blow (MemoryError) up later?
        if function_code is not None and not self.optimise_parameters["simplify expressions"]:
#        if function_code is not None:
            return function_code

        # Get auxiliary variables to generate function
        component, local_comp, local_offset, ffc_element, quad_element, \
        transformation, multiindices = self._get_auxiliary_variables(o, components, derivatives)


        # Create code for function and add empty tuple to cache dict.
        function_code = {(): self.create_function(o, derivatives, component,
                              local_comp, local_offset, ffc_element, quad_element,
                              transformation, multiindices)}

        self.function_cache[(o, components, derivatives, self.restriction)] = function_code

        return function_code

    def constant(self, o, *operands):
        #print("\n\nVisiting Constant: " + repr(o))

        # Safety checks.
        ffc_assert(not operands, "Didn't expect any operands for Constant: " + repr(operands))
        ffc_assert(len(self.component()) == 0, "Constant does not expect component indices: " + repr(self._components))
        ffc_assert(o.shape() == (), "Constant should not have a value shape: " + repr(o.shape()))

        # Component default is 0
        component = 0

        # Handle restriction.
        if self.restriction == "-":
            component += 1

        # Let child class create constant symbol
        coefficient = format["coefficient"](o.count(), component)
        return self._create_symbol(coefficient, CONST)

    def vector_constant(self, o, *operands):
        #print("\n\nVisiting VectorConstant: " + repr(o))

        # Get the component
        components = self.component()

        # Safety checks.
        ffc_assert(not operands, "Didn't expect any operands for VectorConstant: " + repr(operands))
        ffc_assert(len(components) == 1, "VectorConstant expects 1 component index: " + repr(components))

        # We get one component.
        component = components[0]

        # Handle restriction.
        if self.restriction == "-":
            component += o.shape()[0]

        # Let child class create constant symbol
        coefficient = format["coefficient"](o.count(), component)
        return self._create_symbol(coefficient, CONST)

    def tensor_constant(self, o, *operands):
        #print("\n\nVisiting TensorConstant: " + repr(o))

        # Get the components
        components = self.component()

        # Safety checks.
        ffc_assert(not operands, "Didn't expect any operands for TensorConstant: " + repr(operands))
        ffc_assert(len(components) == len(o.shape()), \
                   "The number of components '%s' must be equal to the number of shapes '%s' for TensorConstant." % (repr(components), repr(o.shape())))

        # Let the UFL element handle the component map.
        component = o.element()._sub_element_mapping[components]

        # Handle restriction (offset by value shape).
        if self.restriction == "-":
            component += product(o.shape())

        # Let child class create constant symbol
        coefficient = format["coefficient"](o.count(), component)
        return self._create_symbol(coefficient, CONST)

    # -------------------------------------------------------------------------
    # Indexed (indexed.py).
    # -------------------------------------------------------------------------
    def indexed(self, o):
        #print("\n\nVisiting Indexed:" + repr(o))

        # Get indexed expression and index, map index to current value and update components
        indexed_expr, index = o.operands()
        self._components.push(self.visit(index))

        # Visit expression subtrees and generate code.
        code = self.visit(indexed_expr)

        # Remove component again
        self._components.pop()

        return code

    # -------------------------------------------------------------------------
    # MultiIndex (indexing.py).
    # -------------------------------------------------------------------------
    def multi_index(self, o):
        #print("\n\nVisiting MultiIndex:" + repr(o))

        # Loop all indices in MultiIndex and get current values
        subcomp = []
        for i in o:
            if isinstance(i, FixedIndex):
                subcomp.append(i._value)
            elif isinstance(i, Index):
                subcomp.append(self._index2value[i])

        return tuple(subcomp)

    # -------------------------------------------------------------------------
    # IndexSum (indexsum.py).
    # -------------------------------------------------------------------------
    def index_sum(self, o):
        #print("\n\nVisiting IndexSum: " + str(tree_format(o)))

        # Get expression and index that we're summing over
        summand, multiindex = o.operands()
        index, = multiindex

        # Loop index range, update index/value dict and generate code
        ops = []
        for i in range(o.dimension()):
            self._index2value.push(index, i)
            ops.append(self.visit(summand))
            self._index2value.pop()

        # Call sum to generate summation
        code = self.sum(o, *ops)

        return code

    # -------------------------------------------------------------------------
    # MathFunctions (mathfunctions.py).
    # -------------------------------------------------------------------------
    def sqrt(self, o, *operands):
        #print("\n\nVisiting Sqrt: " + repr(o) + "with operands: " + "\n".join(map(repr,operands)))
        return self._math_function(operands, format["sqrt"])

    def exp(self, o, *operands):
        #print("\n\nVisiting Exp: " + repr(o) + "with operands: " + "\n".join(map(repr,operands)))
        return self._math_function(operands, format["exp"])

    def ln(self, o, *operands):
        #print("\n\nVisiting Ln: " + repr(o) + "with operands: " + "\n".join(map(repr,operands)))
        return self._math_function(operands, format["ln"])

    def cos(self, o, *operands):
        #print("\n\nVisiting Cos: " + repr(o) + "with operands: " + "\n".join(map(repr,operands)))
        return self._math_function(operands, format["cos"])

    def sin(self, o, *operands):
        #print("\n\nVisiting Sin: " + repr(o) + "with operands: " + "\n".join(map(repr,operands)))
        return self._math_function(operands, format["sin"])

    # -------------------------------------------------------------------------
    # PositiveRestricted and NegativeRestricted (restriction.py).
    # -------------------------------------------------------------------------
    def positive_restricted(self, o):
        #print("\n\nVisiting PositiveRestricted: " + repr(o))

        # Just get the first operand, there should only be one.
        restricted_expr = o.operands()
        ffc_assert(len(restricted_expr) == 1, "Only expected one operand for restriction: " + repr(restricted_expr))
        ffc_assert(self.restriction is None, "Expression is restricted twice: " + repr(restricted_expr))

        # Set restriction, visit operand and reset restriction
        self.restriction = "+"
        code = self.visit(restricted_expr[0])
        self.restriction = None

        return code

    def negative_restricted(self, o):
        #print("\n\nVisiting NegativeRestricted: " + repr(o))

        # Just get the first operand, there should only be one.
        restricted_expr = o.operands()
        ffc_assert(len(restricted_expr) == 1, "Only expected one operand for restriction: " + repr(restricted_expr))
        ffc_assert(self.restriction is None, "Expression is restricted twice: " + repr(restricted_expr))

        # Set restriction, visit operand and reset restriction
        self.restriction = "-"
        code = self.visit(restricted_expr[0])
        self.restriction = None

        return code

    # -------------------------------------------------------------------------
    # ComponentTensor (tensors.py).
    # -------------------------------------------------------------------------
    def component_tensor(self, o):
        #print("\n\nVisiting ComponentTensor:\n" + str(tree_format(o)))

        # Get expression and indices
        component_expr, indices = o.operands()

        # Get current component(s)
        components = self.component()

        ffc_assert(len(components) == len(indices), \
                   "The number of known components must be equal to the number of components of the ComponentTensor for this to work.")

        # Update the index dict (map index values of current known indices to
        # those of the component tensor)
        for i, v in izip(indices._indices, components):
            self._index2value.push(i, v)

        # Push an empty component tuple
        self._components.push(())

        # Visit expression subtrees and generate code.
        code = self.visit(component_expr)

        # Remove the index map from the StackDict
        for i in range(len(components)):
            self._index2value.pop()

        # Remove the empty component tuple
        self._components.pop()

        return code

    def list_tensor(self, o):
        #print("\n\nVisiting ListTensor: " + repr(o))

        # Get the component
        component = self.component()

        # Extract first and the rest of the components
        c0, c1 = component[0], component[1:]

        # Get first operand
        op = o.operands()[c0]

        # Evaluate subtensor with this subcomponent
        self._components.push(c1)
        code = self.visit(op)
        self._components.pop()

        return code

    # -------------------------------------------------------------------------
    # Variable (variable.py).
    # -------------------------------------------------------------------------
    def variable(self, o):
        #print("\n\nVisiting Variable: " + repr(o))
        # Just get the expression associated with the variable
        return self.visit(o.expression())

    # -------------------------------------------------------------------------
    # Generate terms for representation.
    # -------------------------------------------------------------------------
    def generate_terms(self, integrand):
        "Generate terms for code generation."
        # Get terms.
        terms = self.visit(integrand)

        f_nzc = format["nonzero columns"](0).split("0")[0]

        # Loop code and add weight and scale factor.
        new_terms = {}
        for key, val in terms.items():
            # If value was zero continue.
            if val is None:
                continue

            # Create data and insert.
            data = self._create_entry_data(val)
            used_nzcs = set([int(k[1].split(f_nzc)[1].split("[")[0]) for k in key if f_nzc in k[1]])
            data.append(used_nzcs)
            new_terms[key] = data
        return new_terms

    # -------------------------------------------------------------------------
    # Helper functions for transformation of UFL objects in base class
    # -------------------------------------------------------------------------
    def _create_symbol(self, symbol, domain):
        error("This function should be implemented by the child class.")

    def _create_product(self, symbols):
        error("This function should be implemented by the child class.")

    def _format_scalar_value(self, value):
        error("This function should be implemented by the child class.")

    def _math_function(self, operands, format_function):
        error("This function should be implemented by the child class.")

    def _get_auxiliary_variables(self, ufl_function, component, derivatives):
        "Helper function for both Coefficient and Argument."

        # Get local component (in case we have mixed elements).
        local_comp, local_elem = ufl_function.element().extract_component(component)

        # Check that we don't take derivatives of QuadratureElements.
        quad_element = local_elem.family() == "Quadrature"
        ffc_assert(not (derivatives and quad_element), \
                   "Derivatives of Quadrature elements are not supported: " + repr(ufl_function))

        # Create FFC element.
        ffc_element = create_element(ufl_function.element())

        # Get relevant sub element and mapping.
        sub_element = create_element(ufl_function.element().extract_component(tuple(component))[1])
        # Assuming that mappings for all basisfunctions are equal (they should be).
        transformation = sub_element.mapping()[0]

        # Handle tensor elements.
        if len(local_comp) > 1:
            local_comp = local_elem._sub_element_mapping[local_comp]
        elif local_comp:
            local_comp = local_comp[0]
        else:
            local_comp = 0

        # Map component using UFL map
        # NOTE: We rely implicitly on a similar ordering of sub elements of mixed elements in UFL and FFC.
        if len(component) > 1:
            component = ufl_function.element()._sub_element_mapping[tuple(component)]
        elif component:
            component = component[0]

        # Compute the local offset (needed for non-affine mappings).
        local_offset = 0
        if component:
            local_offset = component - local_comp

        # Set geo_dim.
        # TODO: All terms REALLY have to be defined on cell with the same
        # geometrical dimension so only do this once and exclude the check?
        geo_dim = ufl_function.element().cell().geometric_dimension()
        if self.geo_dim:
            ffc_assert(geo_dim == self.geo_dim, \
                       "All terms must be defined on cells with the same geometrical dimension.")
        else:
            self.geo_dim = geo_dim

        # Generate FFC multi index for derivatives.
        multiindices = FFCMultiIndex([range(geo_dim)]*len(derivatives)).indices

        return (component, local_comp, local_offset, ffc_element, quad_element, transformation, multiindices)

    def _create_mapping_basis(self, component, deriv, ufl_argument, ffc_element):
        "Create basis name and mapping from given basis_info."

        # Get string for integration points.
        f_ip = format["integration points"]
        generate_psi_name = format["psi name"]

        # Only support test and trial functions.
        # TODO: Verify that test and trial functions will ALWAYS be rearranged to 0 and 1.
        indices = {-2: format["first free index"],
                   -1: format["second free index"],
                    0: format["first free index"],
                    1: format["second free index"]}

        # Check that we have a basis function.
        ffc_assert(ufl_argument.count() in indices, \
                   "Currently, Argument index must be either -2, -1, 0 or 1: " + repr(ufl_argument))

        # Handle restriction through facet.
        facet = {"+": self.facet0, "-": self.facet1, None: self.facet0}[self.restriction]

        # Get element counter and loop index.
        element_counter = self.element_map[self.points][ufl_argument.element()]
        loop_index = indices[ufl_argument.count()]

        # Create basis access, we never need to map the entry in the basis table
        # since we will either loop the entire space dimension or the non-zeros.
        if self.points == 1:
            f_ip = "0"
        basis_access = format["component"]("", [f_ip, loop_index])

        # Offset element space dimension in case of negative restriction,
        # need to use the complete element for offset in case of mixed element.
        space_dim = ffc_element.space_dimension()
        offset = {"+": "", "-": str(space_dim), None: ""}[self.restriction]

        # If we have a restricted function multiply space_dim by two.
        if self.restriction == "+" or self.restriction == "-":
            space_dim *= 2

        # Generate psi name and map to correct values.
        name = generate_psi_name(element_counter, facet, component, deriv)
        name, non_zeros, zeros, ones = self.name_map[name]
        loop_index_range = shape(self.unique_tables[name])[1]

        basis = ""
        # Ignore zeros if applicable
        if zeros and (self.optimise_parameters["ignore zero tables"] or self.optimise_parameters["remove zero terms"]):
            basis = self._format_scalar_value(None)[()]
        # If the loop index range is one we can look up the first component
        # in the psi array. If we only have ones we don't need the basis.
        elif self.optimise_parameters["ignore ones"] and loop_index_range == 1 and ones:
            loop_index = "0"
            basis = self._format_scalar_value(1.0)[()]
        else:
            # Add basis name to the psi tables map for later use.
            basis = self._create_symbol(name + basis_access, BASIS)[()]
            self.psi_tables_map[basis] = name

        # Create the correct mapping of the basis function into the local element tensor.
        basis_map = loop_index
        if non_zeros and basis_map == "0":
            basis_map = str(non_zeros[1][0])
        elif non_zeros:
            basis_map = format["component"](format["nonzero columns"](non_zeros[0]), basis_map)
        if offset:
            basis_map = format["grouping"](format["add"]([basis_map, offset]))

        # Try to evaluate basis map ("3 + 2" --> "5").
        try:
            basis_map = str(eval(basis_map))
        except:
            pass

        # Create mapping (index, map, loop_range, space_dim).
        # Example dx and ds: (0, j, 3, 3)
        # Example dS: (0, (j + 3), 3, 6), 6=2*space_dim
        # Example dS optimised: (0, (nz2[j] + 3), 2, 6), 6=2*space_dim
        mapping = ((ufl_argument.count(), basis_map, loop_index_range, space_dim),)

        return (mapping, basis)

    def _create_function_name(self, component, deriv, quad_element, ufl_function, ffc_element):

        # Get string for integration points.
        f_ip = format["integration points"]
        generate_psi_name = format["psi name"]

        # Pick first free index of secondary type
        # (could use primary indices, but it's better to avoid confusion).
        loop_index = format["free indices"][0]

        # Create basis access, we never need to map the entry in the basis
        # table since we will either loop the entire space dimension or the
        # non-zeros.
        if self.points == 1:
            f_ip = "0"
        basis_access = format["component"]("", [f_ip, loop_index])

        # Handle restriction through facet.
        facet = {"+": self.facet0, "-": self.facet1, None: self.facet0}[self.restriction]

        # Get the element counter.
        element_counter = self.element_map[self.points][ufl_function.element()]

        # Offset by element space dimension in case of negative restriction.
        offset = {"+": "", "-": str(ffc_element.space_dimension()), None: ""}[self.restriction]

        # Create basis name and map to correct basis and get info.
        psi_name = generate_psi_name(element_counter, facet, component, deriv)
        psi_name, non_zeros, zeros, ones = self.name_map[psi_name]

        # If all basis are zero we just return None.
        if zeros and self.optimise_parameters["ignore zero tables"]:
            return self._format_scalar_value(None)[()]

        # Get the index range of the loop index.
        loop_index_range = shape(self.unique_tables[psi_name])[1]

        # Set default coefficient access.
        coefficient_access = loop_index

        # If the loop index range is one we can look up the first component
        # in the coefficient array. If we only have ones we don't need the basis.
        basis_name = psi_name
        if self.optimise_parameters["ignore ones"] and loop_index_range == 1 and ones:
            coefficient_access = "0"
            basis_name = ""
        elif not quad_element:
            # Add basis name to set of used tables and add matrix access.
            # TODO: We should first add this table if the function is used later
            # in the expressions. If some term is multiplied by zero and it falls
            # away there is no need to compute the function value
            self.used_psi_tables.add(psi_name)
            basis_name += basis_access

        # If we have a quadrature element we can use the ip number to look
        # up the value directly. Need to add offset in case of components.
        if quad_element:
            quad_offset = 0
            if component:
                # FIXME: Should we add a member function elements() to FiniteElement?
                if isinstance(ffc_element, MixedElement):
                    for i in range(component):
                        quad_offset += ffc_element.elements()[i].space_dimension()
                elif component != 1:
                    error("Can't handle components different from 1 if we don't have a MixedElement.")
                else:
                    quad_offset += ffc_element.space_dimension()
            if quad_offset:
                coefficient_access = format["add"]([f_ip, str(quad_offset)])
            else:
                coefficient_access = f_ip

        # If we have non zero column mapping but only one value just pick it.
        used_nzcs = set()
        if non_zeros and coefficient_access == "0":
            coefficient_access = str(non_zeros[1][0])
        elif non_zeros and not quad_element:
            used_nzcs.add(non_zeros[0])
            coefficient_access = format["component"](format["nonzero columns"](non_zeros[0]), coefficient_access)
        if offset:
            coefficient_access = format["add"]([coefficient_access, offset])

        # Try to evaluate coefficient access ("3 + 2" --> "5").
        ACCESS = IP
        try:
            coefficient_access = str(eval(coefficient_access))
            ACCESS = GEO
        except:
            pass

        coefficient = format["coefficient"](str(ufl_function.count()), coefficient_access)
        function_expr = self._create_symbol(coefficient, ACCESS)[()]
        if basis_name:
            function_expr = self._create_product([self._create_symbol(basis_name, ACCESS)[()], self._create_symbol(coefficient, ACCESS)[()]])

        # If we have a quadrature element (or if basis was deleted) we don't need the basis.
        if quad_element or not basis_name:
            function_name = self._create_symbol(coefficient, ACCESS)[()]
        else:
            # Check if the expression to compute the function value is already in
            # the dictionary of used function. If not, generate a new name and add.
            function_name = self._create_symbol(format["function value"](self.function_count), ACCESS)[()]
            if not function_expr in self.functions:
                function_name = self._create_symbol(format["function value"](self.function_count), ACCESS)[()]
                self.functions[function_expr] =\
                    (self.function_count, loop_index_range, self._count_operations(function_expr), psi_name, used_nzcs)
                # Increase count.
                self.function_count += 1
            else:
                function_count, index_r, ops, psi_name, used_nzcs = self.functions[function_expr]
                function_name = self._create_symbol(format["function value"](function_count), ACCESS)[()]
                # Check just to make sure.
                ffc_assert(index_r == loop_index_range, "Index ranges does not match." + repr(index_r) + repr(loop_index_range))
        return function_name

    # -------------------------------------------------------------------------
    # Helper functions for code_generation()
    # -------------------------------------------------------------------------
    def _count_operations(self, expression):
        error("This function should be implemented by the child class.")

    def _create_entry_data(self, val):
        error("This function should be implemented by the child class.")
