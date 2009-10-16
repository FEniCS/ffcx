"QuadratureTransformerBase, a common class for quadrature transformers to translate UFL expressions."

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2009-10-13 -- 2009-10-13"
__copyright__ = "Copyright (C) 2009 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Python modules
from itertools import izip

# UFL Classes.
from ufl.classes import MultiIndex
from ufl.classes import FixedIndex
from ufl.classes import Index
from ufl.common import StackDict
from ufl.common import Stack

# UFL Algorithms.
from ufl.algorithms.transformations import Transformer
from ufl.algorithms.transformations import ReuseTransformer

# FFC common modules.
from ffc.common.log import debug, error

# Utility and optimisation functions for quadraturegenerator.
from quadraturegenerator_utils import create_psi_tables

class QuadratureTransformerBase(Transformer):
#class QuadratureTransformerBase(ReuseTransformer):
    "Transform UFL representation to quadrature code."

    def __init__(self, form_representation, domain_type, optimise_options, format):

        Transformer.__init__(self)

        # Save format, optimise_options, weights and fiat_elements_map.
        self.format = format
        self.optimise_options = optimise_options
        self.quadrature_weights = form_representation.quadrature_weights[domain_type]

        # Create containers and variables.
        self.used_psi_tables = set()
        self.psi_tables_map = {}
        self.used_weights = set()
        self.used_nzcs = set()
        self.geo_consts = {}
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
        self.trans_set = set()
        self.element_map, self.name_map, self.unique_tables =\
              create_psi_tables(form_representation.psi_tables[domain_type],\
                                       self.format["epsilon"], self.optimise_options)

        # Cache.
        self.basis_function_cache = {}
        self.function_cache = {}

    def update_facets(self, facet0, facet1):
        self.facet0 = facet0
        self.facet1 = facet1
        # Reset functions and count everytime we generate a new case of facets.
        self.functions = {}
        self.function_count = 0

    def update_points(self, points):
        self.points = points
        # Reset functions everytime we move to a new quadrature loop
        # But not the functions count.
        self.functions = {}

    def reset(self):
        # Reset containers.
        self.used_psi_tables = set()
        self.psi_tables_map = {}
        self.used_weights = set()
        self.used_nzcs = set()
        self.geo_consts = {}
        self.trans_set = set()
        self.functions = {}
        self.function_count = 0
        self.geo_dim = 0
        self.points = 0
        self.facet0 = None
        self.facet1 = None
        if self._components:
            error("This list is supposed to be empty.")
        # It should be zero but clear just to be sure.
        self._components = Stack()
        self._index2value = StackDict()

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
        print "\n\nVisiting basic Expr:", o.__repr__(), "with operands:"
        error("This expression is not handled: ", str(o))

    # Nothing in terminal.py is handled. Can only handle children of these clases.
    def terminal(self, o):
        print "\n\nVisiting basic Terminal:", o.__repr__(), "with operands:"
        error("This terminal is not handled: ", str(o))

    # -------------------------------------------------------------------------
    # Things which should not be here (after expansion etc.) from:
    # algebra.py, differentiation.py, finiteelement.py,
    # form.py, indexing.py, integral.py.
    # -------------------------------------------------------------------------
    def algebra_operator(self, o, *operands):
        print "\n\nVisiting AlgebraOperator: ", o.__repr__()
        error("This type of AlgebraOperator should have been expanded!!" + o.__repr__())

    def derivative(self, o, *operands):
        print "\n\nVisiting Derivative: ", o.__repr__()
        error("All derivatives apart from SpatialDerivative should have been expanded!!")

    def finite_element_base(self, o, *operands):
        print "\n\nVisiting FiniteElementBase: ", o.__repr__()
        error("FiniteElements must be member of a BasisFunction or Function!!")

    def form(self, o, *operands):
        print "\n\nVisiting Form: ", o.__repr__()
        error("The transformer only work on a Form integrand, not the Form itself!!")

    def index_base(self, o):
        print "\n\nVisiting IndexBase: ", o.__repr__()
        error("Indices should not be floating around freely in the integrand!!")

    def integral(self, o):
        print "\n\nVisiting Integral: ", o.__repr__()
        error("Integral should not be present in the integrand!!")

    def measure(self, o):
        print "\n\nVisiting Measure: ", o.__repr__()
        error("Measure should not be present in the integrand!!")

    # -------------------------------------------------------------------------
    # Things which are not supported yet, from:
    # condition.py, constantvalue.py, function.py, geometry.py, mathfunctions.py,
    # restriction.py, tensoralgebra.py, variable.py.
    # -------------------------------------------------------------------------
    def condition(self, o):
        print "\n\nVisiting Condition:", o.__repr__()
        error("Condition is not supported (yet).")

    def conditional(self, o):
        print "\n\nVisiting Condition:", o.__repr__()
        error("Conditional is not supported (yet).")

    def scalar_something(self, o):
        print "\n\nVisiting ConstantSomething:", o.__repr__()
        error("ConstantSomething is not supported (yet).")

    def zero(self, o):
        print "\n\nVisiting Zero:", o.__repr__()
        error("Zero is not supported (yet).")

    def geometric_quantity(self, o):
        print "\n\nVisiting GeometricQuantity:", o.__repr__()
        error("GeometricQuantity is not supported (yet).")

    def spatial_coordinate(self, o):
        print "\n\nVisiting SpatialCoordinate:", o.__repr__()
        error("SpatialCoordinate is not supported (yet).")

    def math_function(self, o):
        print "\n\nVisiting MathFunction:", o.__repr__()
        error("This MathFunction is not supported (yet).")

    def restricted(self, o):
        print "\n\nVisiting Restricted:", o.__repr__()
        error("This type of Restricted is not supported (only positive and negative are supported).")

    # -------------------------------------------------------------------------
    # Things that should be implemented by child classes.
    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------
    # AlgebraOperators (algebra.py).
    # -------------------------------------------------------------------------
    def sum(self, o, *operands):
        print "\n\nVisiting Sum: ", o.__repr__()
        error("This object should be implemented by the child class.")

    def product(self, o, *operands):
        print "\n\nVisiting Product: ", o.__repr__()
        error("This object should be implemented by the child class.")

    def division(self, o, *operands):
        print "\n\nVisiting Division: ", o.__repr__()
        error("This object should be implemented by the child class.")

    def power(self, o):
        print "\n\nVisiting Power: ", o.__repr__()
        error("This object should be implemented by the child class.")

    def abs(self, o, *operands):
        print "\n\nVisiting Abs: ", o.__repr__()
        error("This object should be implemented by the child class.")

    # -------------------------------------------------------------------------
    # Constant values (constantvalue.py).
    # -------------------------------------------------------------------------
    def identity(self, o):
        print "\n\nVisiting Identity: ", o.__repr__()
        error("This object should be implemented by the child class.")

    def scalar_value(self, o, *operands):
        print "\n\nVisiting ScalarValue: ", o.__repr__()
        error("This object should be implemented by the child class.")

    # -------------------------------------------------------------------------
    # Function and Constants (function.py).
    # -------------------------------------------------------------------------
    def constant(self, o, *operands):
        print "\n\nVisiting Constant: ", o.__repr__()
        error("This object should be implemented by the child class.")

    def vector_constant(self, o, *operands):
        print "\n\nVisiting VectorConstant: ", o.__repr__()
        error("This object should be implemented by the child class.")

    def tensor_constant(self, o, *operands):
        print "\n\nVisiting TensorConstant: ", o.__repr__()
        error("This object should be implemented by the child class.")

    # -------------------------------------------------------------------------
    # FacetNormal (geometry.py).
    # -------------------------------------------------------------------------
    def facet_normal(self, o,  *operands):
        print "\n\nVisiting FacetNormal: ", o.__repr__()
        error("This object should be implemented by the child class.")

    # -------------------------------------------------------------------------
    # Things that can be handled by the base class.
    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------
    # BasisFunction (basisfunction.py).
    # -------------------------------------------------------------------------
    def basis_function(self, o, *operands):
        #print("\nVisiting BasisFunction:" + str(o))

        # Just checking that we don't get any operands.
        if operands:
            error("Didn't expect any operands for BasisFunction: " + str(operands))

        # Create aux. info.
        component = self.component()
        derivatives = self.derivatives()

#        print ("BasisFunction: component: " + str(component))
#        print ("BasisFunction: derivatives: " + str(derivatives))

        # Check if basis is already in cache
        basis = self.basis_function_cache.get((o, component, derivatives), None)
        if basis is not None:
            return basis

        # Create mapping and code for basis function.
        basis = self.create_basis_function(o, component, derivatives)
        self.basis_function_cache[(o, component, derivatives)] = basis

#        print "BasisFunction: code: " + str(basis)

        return basis

    # -------------------------------------------------------------------------
    # SpatialDerivative (differentiation.py).
    # -------------------------------------------------------------------------
    def spatial_derivative(self, o):
#        #print("\n\nVisiting SpatialDerivative: " + o.__repr__())

        derivative_expr, index = o.operands()

#        print "\nSpatialDerivative: derivative_expr: ", derivative_expr
#        print "SpatialDerivative: index: ", index

        # Get direction of derivative and check that we only get one return index
        der = self.visit(index)
        if len(der) != 1:
            error("SpatialDerivative: expected only one direction index. " + str(der))

        self._derivatives.append(der[0])
#        print "SpatialDerivative: _directions: ", self._derivatives

        # Visit children to generate the derivative code.
        code = self.visit(derivative_expr)

        self._derivatives.pop()

#        print "SpatialDerivative: code: ", code

        return code

    # -------------------------------------------------------------------------
    # Function (function.py).
    # -------------------------------------------------------------------------
    def function(self, o, *operands):
        #print("\nVisiting Function: " + str(o))

        # Safety check.
        if operands:
            error("Didn't expect any operands for Function: " + str(operands))

        # Create aux. info.
        component = self.component()
        derivatives = self.derivatives()

#        print("component: " + str(component))
#        print("derivatives: " + str(derivatives))

        # Check if function is already in cache
        function_code = self.function_cache.get((o, component, derivatives), None)
        if function_code is not None:
#            print "function cache!"
            return function_code

        # Create code for function and add empty tuple to cache dict.
        function_code = {(): self.create_function(o, component, derivatives)}

        self.function_cache[(o, component, derivatives)] = function_code

        return function_code

    # -------------------------------------------------------------------------
    # Indexed (indexed.py).
    # -------------------------------------------------------------------------
    def indexed(self, o):
#        #print("\n\nVisiting Indexed:" + o.__repr__())

        indexed_expr, index = o.operands()
        self._components.push(self.visit(index))
#        print "\nIndexed:"
#        print "Indexed: indexed_expr: ", indexed_expr
#        print "Indexed: index: ", repr(index)
#        print "Indexed: comps: ", self._components

#        #print("\nwrap, indexed_expr: " + indexed_expr.__repr__())
#        #print("\nwrap, index.__repr__(): " + index.__repr__())

#        # Save copy of components to let parent delete them again.
#        old_components = self._components[:]
#        self._components = []

#        # Loop multi indices and create components.
#        for indx in index:
#            self._components.append(indx)

#        #print("\nnew components: " + str(self._components))

        # Visit expression subtrees and generate code.
        code = self.visit(indexed_expr)

#        #print("\nindexed code: " + str(code))

#        # Set components equal to old components.
#        self._components = old_components[:]
        self._components.pop()

#        #print("\nold components: " + str(self._components))

        return code

    # -------------------------------------------------------------------------
    # MultiIndex (indexing.py).
    # -------------------------------------------------------------------------
    def multi_index(self, o):
#        #print("\n\nVisiting MultiIndex:" + o.__repr__())
#        print "\nMultiIndex: o: ", repr(o)
#        print "\nMultiIndex: o: ", repr(operands)
        subcomp = []
        for i in o:
            if isinstance(i, FixedIndex):
                subcomp.append(i._value)
            elif isinstance(i, Index):
                subcomp.append(self._index2value[i])
#        print "MultiIndex: subcomp: ", tuple(subcomp)
        return tuple(subcomp)

    # -------------------------------------------------------------------------
    # IndexSum (indexsum.py).
    # -------------------------------------------------------------------------
    def index_sum(self, o):
#        #print("\n\nVisiting IndexSum: " + o.__repr__())

        summand, multiindex = o.operands()
        index, = multiindex

#        print "\nIndexSum: summand: ", summand
#        print "IndexSum: multiind: ", repr(multiindex)
#        print "IndexSum: index: ", repr(index)
#        print "IndexSum: o.dim: ", o.dimension()

        ops = []
        for i in range(o.dimension()):
#            print "IndexSum: i: ", i
            self._index2value.push(index, i)
#            print "IndexSum: ind2val: ", self._index2value
            ops.append(self.visit(summand))
            self._index2value.pop()

        code = self.sum(o, *ops)

#        print "IndexSum: code: ", code
        return code

    # -------------------------------------------------------------------------
    # MathFunctions (mathfunctions.py).
    # -------------------------------------------------------------------------
    def sqrt(self, o, *operands):
#        #print("\n\nVisiting Sqrt: " + o.__repr__() + "with operands: " + "\n".join(map(str,operands)))
        # Call common math function.
        return self._math_function(operands, self.format["sqrt"])

    def exp(self, o, *operands):
#        #print("\n\nVisiting Exp: " + o.__repr__() + "with operands: " + "\n".join(map(str,operands)))
        # Call common math function.
        return self._math_function(operands, self.format["exp"])

    def ln(self, o, *operands):
#        #print("\n\nVisiting Ln: " + o.__repr__() + "with operands: " + "\n".join(map(str,operands)))
        # Call common math function.
        return self._math_function(operands, self.format["ln"])

    def cos(self, o, *operands):
#        #print("\n\nVisiting Cos: " + o.__repr__() + "with operands: " + "\n".join(map(str,operands)))
        # Call common math function.
        return self._math_function(operands, self.format["cos"])

    def sin(self, o, *operands):
#        #print("\n\nVisiting Sin: " + o.__repr__() + "with operands: " + "\n".join(map(str,operands)))
        # Call common math function.
        return self._math_function(operands, self.format["sin"])

    # -------------------------------------------------------------------------
    # PositiveRestricted and NegativeRestricted (restriction.py).
    # -------------------------------------------------------------------------
    def positive_restricted(self, o):
#        #print("\n\nVisiting PositiveRestricted: " + o.__repr__())

        # Just get the first operand, there should only be one.
        restricted_expr = o.operands()
        if len(restricted_expr) != 1:
            error("Only expected one operand for restriction: " + str(restricted_expr))
 
        # Visit operand and generate restricted code.
        self.restriction = "+"
        code = self.visit(restricted_expr[0])

        # Reset restriction.
        # TODO: Is this really necessary?
        self.restriction = None

        return code

    def negative_restricted(self, o):
#        #print("\n\nVisiting NegativeRestricted: " + o.__repr__())

        # Just get the first operand, there should only be one.
        restricted_expr = o.operands()
        if len(restricted_expr) != 1:
            error("Only expected one operand for restriction: " + str(restricted_expr))
 
        # Visit operand and generate restricted code.
        self.restriction = "-"
        code = self.visit(restricted_expr[0])

        # Reset restriction.
        # TODO: Is this really necessary?
        self.restriction = None

        return code

    # -------------------------------------------------------------------------
    # ComponentTensor (tensors.py).
    # -------------------------------------------------------------------------
    def component_tensor(self, o):
#        #print("\n\nVisiting ComponentTensor: " + o.__repr__())

        component_expr, indices = o.operands()
#        #print("component_expr: " + str(component_expr))
#        #print("indices: " + repr(indices))

        # Get current component(s)
        comps = self.component()

#        print "\nComponentTensor: component_expr: ", component_expr
#        print "ComponentTensor: indices: ", repr(indices)
#        print "ComponentTensor: comps: ", comps

        if not len(comps) == len(indices):
            error("The number of known components must be equal to the number of components of the ComponentTensor for this to work.")

        # Update the index dict (map index values of current known indices to
        # those of the component tensor)
        for i, v in izip(indices._indices, comps):
            self._index2value.push(i, v)

        # Push an empty component tuple
        self._components.push(())

#        print "ComponentTensor: _components: ", self._components

        # Visit expression subtrees and generate code.
        code = self.visit(component_expr)

        # Remove the index map from the StackDict
        for i in range(len(comps)):
            self._index2value.pop()

        # Remove the empty component tuple
        self._components.pop()

#        print "ComponentTensor: _components (after pop): ", self._components

#        print "ComponentTensor: code: ", code

        return code

    def list_tensor(self, o):
#        #print("\n\nVisiting ListTensor: " + o.__repr__())

#        print "\nListTensor: o: ", repr(o)
        component = self.component()
#        print "ListTensor: component: ", component
#        print "ListTensor: operands: ", o.operands()

#        c = self.component()
        c0, c1 = component[0], component[1:]
#        print "ListTensor: c0 ", c0
#        print "ListTensor: c1 ", c1

        # Get first operand
        op = o.operands()[c0]

        # Evaluate subtensor with this subcomponent
        self._components.push(c1)
        code = self.visit(op)
        self._components.pop()

#        print "ListTensor: code ", code
        return code

#        return r

    # -------------------------------------------------------------------------
    # Variable (variable.py).
    # -------------------------------------------------------------------------
    def variable(self, o):
#        #print("\n\nVisiting Variable: " + o.__repr__())
        return self.visit(o.expression())

