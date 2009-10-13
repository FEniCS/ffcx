"QuadratureTransformerBase, a common class for quadrature transformers to translate UFL expressions."

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2009-10-13 -- 2009-10-13"
__copyright__ = "Copyright (C) 2009 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# UFL Classes.
from ufl.classes import MultiIndex
from ufl.classes import FixedIndex

# UFL Algorithms.
from ufl.algorithms.transformations import Transformer

# FFC common modules.
from ffc.common.log import debug, error

# Utility and optimisation functions for quadraturegenerator.
from quadraturegenerator_utils import create_psi_tables

class QuadratureTransformerBase(Transformer):
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
        self._components = []
        self.trans_set = set()
        self.element_map, self.name_map, self.unique_tables =\
              create_psi_tables(form_representation.psi_tables[domain_type],\
                                       self.format["epsilon"], self.optimise_options)

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
        self._components = []

    def disp(self):
        print "\n\n **** Displaying QuadratureTransformer ****"
        print "\nQuadratureTransformer, element_map:\n", self.element_map
        print "\nQuadratureTransformer, name_map:\n", self.name_map
        print "\nQuadratureTransformer, unique_tables:\n", self.unique_tables
        print "\nQuadratureTransformer, used_psi_tables:\n", self.used_psi_tables
        print "\nQuadratureTransformer, psi_tables_map:\n", self.psi_tables_map
        print "\nQuadratureTransformer, used_weights:\n", self.used_weights
        print "\nQuadratureTransformer, geo_consts:\n", self.geo_consts

    # -------------------------------------------------------------------------
    # Start handling UFL classes.
    # -------------------------------------------------------------------------
    # Nothing in expr.py is handled. Can only handle children of these clases.
    def expr(self, o, *operands):
        print "\n\nVisiting basic Expr:", o.__repr__(), "with operands:"
        print ", ".join(map(str,operands))
        error("This expression is not handled: ", str(o))

    # Nothing in terminal.py is handled. Can only handle children of these clases.
    def terminal(self, o, *operands):
        print "\n\nVisiting basic Terminal:", o.__repr__(), "with operands:"
        print ", ".join(map(str,operands))
        error("This terminal is not handled: ", str(o))

    # -------------------------------------------------------------------------
    # Things which should not be here (after expansion etc.) from:
    # algebra.py, constantvalue.py, differentiation.py, finiteelement.py,
    # form.py, indexing.py, integral.py, tensors.py.
    # -------------------------------------------------------------------------
    def algebra_operator(self, o, *operands):
        debug("\n\nVisiting AlgebraOperator: " + o.__repr__())
        error("This type of AlgebraOperator should have been expanded!!" + o.__repr__())

    def identity(self, o, *operands):
        debug("\n\nVisiting Identity: " + o.__repr__())
        error("Identity should have been expanded!!")

    def derivative(self, o, *operands):
        debug("\n\nVisiting Derivative: " + o.__repr__())
        error("All derivatives apart from SpatialDerivative should have been expanded!!")

    def finite_element_base(self, o, *operands):
        debug("\n\nVisiting FiniteElementBase: " + o.__repr__())
        error("FiniteElements must be member of a BasisFunction or Function!!")

    def form(self, o, *operands):
        debug("\n\nVisiting Form: " + o.__repr__())
        error("The transformer only work on a Form integrand, not the Form itself!!")

    def index_base(self, o):
        debug("\n\nVisiting IndexBase: " + o.__repr__())
        error("Indices should not be floating around freely in the integrand!!")

    def integral(self, o):
        debug("\n\nVisiting Integral: " + o.__repr__())
        error("Integral should not be present in the integrand!!")

    def measure(self, o):
        debug("\n\nVisiting Measure: " + o.__repr__())
        error("Measure should not be present in the integrand!!")

    def list_tensor(self, o):
        debug("\n\nVisiting ListTensor: " + o.__repr__())
        error("ListTensor should have been purged!!")

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

    def variable(self, o):
        print "\n\nVisiting Variable:", o.__repr__()
        error("Variable is not handled yet, should have been removed by strip_variables).")

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
        debug("\n\nVisiting Division: " + o.__repr__() + "with operands: " + "\n".join(map(str,operands)))

    def power(self, o):
        print "\n\nVisiting Power: ", o.__repr__()
        error("This object should be implemented by the child class.")

    def abs(self, o, *operands):
        print "\n\nVisiting Abs: ", o.__repr__()
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
    # Supported classes.
    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------
    # BasisFunction (basisfunction.py).
    # -------------------------------------------------------------------------
    def basis_function(self, o, *operands):
        debug("\n\nVisiting BasisFunction:" + o.__repr__())

        # Just checking that we don't get any operands.
        if operands:
            error("Didn't expect any operands for BasisFunction: " + str(operands))

        # Create aux. info.
        component = []
        derivatives = ()

        # Handle derivatives and components.
        if self._derivatives:
            derivatives = self._derivatives[:]
        if self._components:
            component = [int(c) for c in self._components]
        debug("\ncomponent: " + str(component))
        debug("\nderivatives: " + str(derivatives))

        # Create mapping and code for basis function.
        basis = self.create_basis_function(o, component, derivatives)

        # Reset spatial derivatives.
        # FIXME: (should this be handled by SpatialDerivative).
        self._derivatives = []

        debug("basis code: " + str(basis))

        return basis

    # -------------------------------------------------------------------------
    # SpatialDerivative (differentiation.py).
    # -------------------------------------------------------------------------
    def spatial_derivative(self, o):
        debug("\n\nVisiting SpatialDerivative: " + o.__repr__())

        derivative_expr, index = o.operands()
        debug("derivative_expr: " + str(derivative_expr))
        debug("index: " + index.__repr__())
        if not isinstance(index, MultiIndex) and len(index) == 1 and isinstance(index[0], FixedIndex):
            error("Expecting 1 MultiIndex with 1 FixedIndex: " + str(index))

        # Get the direction that we need the derivative for.
        direction = int(index[0])

        # Append the derivative.
        self._derivatives.append(direction)

        # Visit children to generate the derivative code.
        code = self.visit(derivative_expr)
        debug("Spatial_derivative, code: " + str(code))
        return code

    # -------------------------------------------------------------------------
    # Function (function.py).
    # -------------------------------------------------------------------------
    def function(self, o, *operands):
        debug("\n\nVisiting Function: " + o.__repr__())

        # Safety checks.
        if operands:
            error("Didn't expect any operands for Function: " + str(operands))
        if not all(isinstance(c, FixedIndex) for c in self._components):
            error("Function expects FixedIndex for its components: " + str(self._components))

        # Create aux. info.
        component = []
        derivatives = ()

        # Handle derivatives and components.
        if self._derivatives:
            derivatives = self._derivatives[:]
        if self._components:
            component = [int(c) for c in self._components]

        debug("\ncomponent: " + str(component))
        debug("\nderivatives: " + str(derivatives))

        # Create code for basis function.
        code = self.create_function(o, component, derivatives)

        # Reset spatial derivatives.
        # FIXME: (should this be handled by SpatialDerivative).
        self._derivatives = []

        debug("function code: " + str(code))

        return {(): code}

    # -------------------------------------------------------------------------
    # Indexed (indexed.py).
    # -------------------------------------------------------------------------
    def indexed(self, o):
        debug("\n\nVisiting Indexed:" + o.__repr__())

        indexed_expr, index = o.operands()

        debug("\nwrap, indexed_expr: " + indexed_expr.__repr__())
        debug("\nwrap, index.__repr__(): " + index.__repr__())

        # Save copy of components to let parent delete them again.
        old_components = self._components[:]
        self._components = []

        # Loop multi indices and create components.
        for indx in index:
            self._components.append(indx)

        debug("\nnew components: " + str(self._components))

        # Visit expression subtrees and generate code.
        code = self.visit(indexed_expr)

        debug("\nindexed code: " + str(code))

        # Set components equal to old components.
        self._components = old_components[:]

        debug("\nold components: " + str(self._components))

        return code

    # -------------------------------------------------------------------------
    # MathFunctions (mathfunctions.py).
    # -------------------------------------------------------------------------
    def sqrt(self, o, *operands):
        debug("\n\nVisiting Sqrt: " + o.__repr__() + "with operands: " + "\n".join(map(str,operands)))
        # Call common math function.
        return self._math_function(operands, self.format["sqrt"])

    def exp(self, o, *operands):
        debug("\n\nVisiting Exp: " + o.__repr__() + "with operands: " + "\n".join(map(str,operands)))
        # Call common math function.
        return self._math_function(operands, self.format["exp"])

    def ln(self, o, *operands):
        debug("\n\nVisiting Ln: " + o.__repr__() + "with operands: " + "\n".join(map(str,operands)))
        # Call common math function.
        return self._math_function(operands, self.format["ln"])

    def cos(self, o, *operands):
        debug("\n\nVisiting Cos: " + o.__repr__() + "with operands: " + "\n".join(map(str,operands)))
        # Call common math function.
        return self._math_function(operands, self.format["cos"])

    def sin(self, o, *operands):
        debug("\n\nVisiting Sin: " + o.__repr__() + "with operands: " + "\n".join(map(str,operands)))
        # Call common math function.
        return self._math_function(operands, self.format["sin"])

    # -------------------------------------------------------------------------
    # PositiveRestricted and NegativeRestricted (restriction.py).
    # -------------------------------------------------------------------------
    def positive_restricted(self, o):
        debug("\n\nVisiting PositiveRestricted: " + o.__repr__())

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

        debug("code: " + str(code))

        return code

    def negative_restricted(self, o):
        debug("\n\nVisiting NegativeRestricted: " + o.__repr__())

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

        debug("code: " + str(code))

        return code

    # -------------------------------------------------------------------------
    # ComponentTensor (tensors.py).
    # -------------------------------------------------------------------------
    def component_tensor(self, o):
        debug("\n\nVisiting ComponentTensor: " + o.__repr__())

        component_expr, index = o.operands()
        debug("component_expr: " + str(component_expr))
        debug("index: " + index.__repr__())

        if not len(self._components) == len(index):
            error("The number of known components must be equal to the number of components of the ComponentTensor for this to work.")

        # Save copy of components to let parent delete them again.
        old_components = self._components[:]
        self._components = []

        debug("\nnew components: " + str(self._components))

        # Visit expression subtrees and generate code.
        code = self.visit(component_expr)

        debug("code: " + str(code))

        # Set components equal to old components.
        self._components = old_components[:]

        debug("\nold components: " + str(self._components))

        return code

    # -------------------------------------------------------------------------
    # Variable (variable.py).
    # -------------------------------------------------------------------------
#    def variable(self, o):
#        debug("\n\nVisiting Variable: " + o.__repr__())
#        return self.visit(o.expression())

