"QuadratureTransformerBase, a common class for quadrature transformers to translate UFL expressions."

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2009-10-13 -- 2009-10-19"
__copyright__ = "Copyright (C) 2009 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Python modules
from itertools import izip
import time

# UFL Classes.
from ufl.classes import MultiIndex
from ufl.classes import FixedIndex
from ufl.classes import Index
from ufl.common import StackDict
from ufl.common import Stack

# UFL Algorithms.
from ufl.algorithms import propagate_restrictions
from ufl.algorithms.transformations import Transformer
from ufl.algorithms.printing import tree_format

# FFC common modules.
from ffc.common.log import debug, error, info

# FFC compiler modules.
from ffc.compiler.tensor.multiindex import MultiIndex as FFCMultiIndex

# FFC fem modules.
from ffc.fem.createelement import create_element

# Utility and optimisation functions for quadraturegenerator.
from quadraturegenerator_utils import create_psi_tables
from quadraturegenerator_utils import generate_loop
from symbolics import generate_aux_constants

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

        # Reset cache
        self.basis_function_cache = {}
        self.function_cache = {}

    def update_points(self, points):
        self.points = points
        # Reset functions everytime we move to a new quadrature loop
        # But not the functions count.
        self.functions = {}

        # Reset cache
        self.basis_function_cache = {}
        self.function_cache = {}

    def reset(self):
        # Reset containers.
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
        if self._components:
            error("This list is supposed to be empty.")
        # It should be zero but clear just to be sure.
        self._components = Stack()
        self._index2value = StackDict()

        # Reset cache
        self.basis_function_cache = {}
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
        print "\n\nVisiting basic Expr:", o.__repr__(), "with operands:"
        error("This expression is not handled: ", str(o))

    # Nothing in terminal.py is handled. Can only handle children of these clases.
    def terminal(self, o):
        print "\n\nVisiting basic Terminal:", o.__repr__(), "with operands:"
        error("This terminal is not handled: ", str(o))

    # -------------------------------------------------------------------------
    # Things which should not be here (after expansion etc.) from:
    # algebra.py, differentiation.py, finiteelement.py,
    # form.py, geometry.py, indexing.py, integral.py, tensoralgebra.py, variable.py.
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

    def space(self, o):
        print "\n\nVisiting Space: ", o.__repr__()
        error("A Space should not be present in the integrand.")

    def cell(self, o):
        print "\n\nVisiting Cell: ", o.__repr__()
        error("A Cell should not be present in the integrand.")

    def index_base(self, o):
        print "\n\nVisiting IndexBase: ", o.__repr__()
        error("Indices should not be floating around freely in the integrand!!")

    def integral(self, o):
        print "\n\nVisiting Integral: ", o.__repr__()
        error("Integral should not be present in the integrand!!")

    def measure(self, o):
        print "\n\nVisiting Measure: ", o.__repr__()
        error("Measure should not be present in the integrand!!")

    def compound_tensor_operator(self, o):
        print "\n\nVisiting CompoundTensorOperator: ", o.__repr__()
        error("CompoundTensorOperator should have been expanded.")

    def label(self, o):
        print "\n\nVisiting Label: ", o.__repr__()
        error("What is a Lable doing in the integrand?")

    # -------------------------------------------------------------------------
    # Things which are not supported yet, from:
    # condition.py, constantvalue.py, function.py, geometry.py, lifting.py,
    # mathfunctions.py, restriction.py
    # -------------------------------------------------------------------------
    def condition(self, o):
        print "\n\nVisiting Condition:", o.__repr__()
        error("Condition is not supported (yet).")

    def conditional(self, o):
        print "\n\nVisiting Condition:", o.__repr__()
        error("Conditional is not supported (yet).")

    def constant_value(self, o):
        print "\n\nVisiting ConstantValue:", o.__repr__()
        error("This type of ConstantValue is not supported (yet).")

    def index_annotated(self, o):
        print "\n\nVisiting IndexAnnotated:", o.__repr__()
        error("Only child classes of IndexAnnotated is supported.")

    def constant_base(self, o):
        print "\n\nVisiting ConstantBase:", o.__repr__()
        error("This type of ConstantBase is not supported (yet).")

    def geometric_quantity(self, o):
        print "\n\nVisiting GeometricQuantity:", o.__repr__()
        error("This type of GeometricQuantity is not supported (yet).")

    def spatial_coordinate(self, o):
        print "\n\nVisiting SpatialCoordinate:", o.__repr__()
        error("SpatialCoordinate is not supported (yet).")

    def lifting_result(self, o):
        print "\n\nVisiting LiftingResult:", o.__repr__()
        error("LiftingResult (and children) is not supported (yet).")

    def terminal_operator(self, o):
        print "\n\nVisiting TerminalOperator:", o.__repr__()
        error("TerminalOperator (LiftingOperator and LiftingFunction) is not supported (yet).")

    def math_function(self, o):
        print "\n\nVisiting MathFunction:", o.__repr__()
        error("This MathFunction is not supported (yet).")

    def restricted(self, o):
        print "\n\nVisiting Restricted:", o.__repr__()
        error("This type of Restricted is not supported (only positive and negative are currently supported).")

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
    # FacetNormal (geometry.py).
    # -------------------------------------------------------------------------
    def facet_normal(self, o,  *operands):
        print "\n\nVisiting FacetNormal: ", o.__repr__()
        error("This object should be implemented by the child class.")

    # -------------------------------------------------------------------------
    # Common auxiliary functions.
    # -------------------------------------------------------------------------
    def get_auxiliary_variables(self, ufl_function, component, derivatives):
        "Helper function for both Function and BasisFunction."

        # Get local component (in case we have mixed elements).
        local_comp, local_elem = ufl_function.element().extract_component(component)

        # Check that we don't take derivatives of QuadratureElements.
        quad_element = local_elem.family() == "Quadrature"
        if derivatives and quad_element:
            error("Derivatives of Quadrature elements are not supported: " + str(ufl_function))

        # Handle tensor elements.
        if len(local_comp) > 1:
            local_comp = local_elem._sub_element_mapping[local_comp]
        elif local_comp:
            local_comp = local_comp[0]
        else:
            local_comp = 0

        # Map component
        if len(component) > 1:
            component = ufl_function.element()._sub_element_mapping[tuple(component)]
        elif component:
            component = component[0]

        # Compute the local offset (needed for non-affine mappings).
        local_offset = 0
        if component:
            local_offset = component - local_comp

        # Create FFC element and get transformation.
        ffc_element = create_element(ufl_function.element())
        transformation = ffc_element.component_element(component)[0].mapping()

        # Set geo_dim.
        # TODO: All terms REALLY have to be defined on cell with the same
        # geometrical dimension so only do this once and exclude the check?
        geo_dim = ufl_function.element().cell().geometric_dimension()
        if self.geo_dim:
            if geo_dim != self.geo_dim:
                error("All terms must be defined on cells with the same geometrical dimension.")
        else:
            self.geo_dim = geo_dim

        # Generate FFC multi index for derivatives.
        multiindices = FFCMultiIndex([range(geo_dim)]*len(derivatives)).indices

        return (component, local_comp, local_offset, ffc_element, quad_element, transformation, multiindices)

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
        components = self.component()
        derivatives = self.derivatives()

        #print ("BasisFunction: components: " + str(components))
        #print ("BasisFunction: derivatives: " + str(derivatives))

        # Check if basis is already in cache
        basis = self.basis_function_cache.get((o, components, derivatives, self.restriction), None)
        # FIXME: Why does using a code dict from cache make the expression manipulations blow (MemoryError) up later?
        if basis is not None and not self.optimise_options["simplify expressions"]:
#        if basis is not None:
            return basis

        # Get auxiliary variables to generate basis
        component, local_comp, local_offset, ffc_element, quad_element, \
        transformation, multiindices = self.get_auxiliary_variables(o, components, derivatives)

        # Create mapping and code for basis function and add to dict.
        basis = self.create_basis_function(o, derivatives, component, local_comp,
                  local_offset, ffc_element, transformation, multiindices)

        self.basis_function_cache[(o, components, derivatives, self.restriction)] = basis

        return basis

    # -------------------------------------------------------------------------
    # Constant values (constantvalue.py).
    # -------------------------------------------------------------------------
    def identity(self, o):
        #print "\n\nVisiting Identity: ", o.__repr__()

        # Get components
        components = self.component()

        # Safety checks.
        if o.operands():
            error("Didn't expect any operands for Identity: " + str(o.operands()))
        elif len(components) != 2:
            error("Identity expect exactly two component indices: " + str(components))

        # Only return a value if i==j
        if components[0] == components[1]:
            return self.format_scalar_value(1.0)
        return self.format_scalar_value(None)

    def scalar_value(self, o, *operands):
        "ScalarValue covers IntValue and FloatValue"
        #print "\n\nVisiting ScalarValue: ", o.__repr__()

        # FIXME: Might be needed because it can be IndexAnnotated?
        if operands:
            error("Did not expect any operands for ScalarValue: " + str((o, operands)))

        return self.format_scalar_value(o.value())

    def zero(self, o, *operands):
        #print "\n\nVisiting Zero:", o.__repr__()
        # FIXME: Might be needed because it can be IndexAnnotated?
        if operands:
            error("Did not expect any operands for Zero: " + str((o, operands)))
        return self.format_scalar_value(None)

    # -------------------------------------------------------------------------
    # SpatialDerivative (differentiation.py).
    # -------------------------------------------------------------------------
    def spatial_derivative(self, o):
        #print("\n\nVisiting SpatialDerivative: " + o.__repr__())

        # Get expression and index
        derivative_expr, index = o.operands()

        # Get direction of derivative and check that we only get one return index
        der = self.visit(index)
        if len(der) != 1:
            error("SpatialDerivative: expected only one direction index. " + str(der))

        # Add direction to list of derivatives
        self._derivatives.append(der[0])

        # Visit children to generate the derivative code.
        code = self.visit(derivative_expr)

        # Remove the direction from list of derivatives
        self._derivatives.pop()
        return code

    # -------------------------------------------------------------------------
    # Function and Constants (function.py).
    # -------------------------------------------------------------------------
    def function(self, o, *operands):
        #print("\nVisiting Function: " + repr(o))

        # Safety check.
        if operands:
            error("Didn't expect any operands for Function: " + str(operands))

        # Create aux. info.
        components = self.component()
        derivatives = self.derivatives()

        #print("components: " + str(components))
        #print("derivatives: " + str(derivatives))

        # Check if function is already in cache
        function_code = self.function_cache.get((o, components, derivatives, self.restriction), None)
        # FIXME: Why does using a code dict from cache make the expression manipulations blow (MemoryError) up later?
#        if function_code is not None:
        if function_code is not None and not self.optimise_options["simplify expressions"]:
            return function_code

        # Get auxiliary variables to generate function
        component, local_comp, local_offset, ffc_element, quad_element, \
        transformation, multiindices = self.get_auxiliary_variables(o, components, derivatives)


        # Create code for function and add empty tuple to cache dict.
        function_code = {(): self.create_function(o, derivatives, component,
                              local_comp, local_offset, ffc_element, quad_element,
                              transformation, multiindices)}

        self.function_cache[(o, components, derivatives, self.restriction)] = function_code

        return function_code

    def constant(self, o, *operands):
        #print("\n\nVisiting Constant: " + o.__repr__())

        # Safety checks.
        if operands:
            error("Didn't expect any operands for Constant: " + str(operands))
        elif len(self.component()) > 0:
            error("Constant does not expect component indices: " + str(self._components))
        elif o.shape() != ():
            error("Constant should not have a value shape: " + str(o.shape()))

        # Component default is 0
        component = 0

        # Handle restriction.
        if self.restriction == "-":
            component += 1

        # Let child class handle the constant coefficient
        return self.create_constant_coefficient(o.count(), component)

    def vector_constant(self, o, *operands):
        #print("\n\nVisiting VectorConstant: " + o.__repr__())

        # Get the component
        components = self.component()

        # Safety checks.
        if operands:
            error("Didn't expect any operands for VectorConstant: " + str(operands))
        elif len(components) != 1:
            error("VectorConstant expects 1 component index: " + str(components))

        # We get one component.
        component = components[0]

        # Handle restriction.
        if self.restriction == "-":
            component += o.shape()[0]

        # Let child class handle the constant coefficient
        return self.create_constant_coefficient(o.count(), component)

    def tensor_constant(self, o, *operands):
        #print("\n\nVisiting TensorConstant: " + o.__repr__())

        # Get the components
        components = self.component()

        # Safety checks.
        if operands:
            error("Didn't expect any operands for TensorConstant: " + str(operands))
        elif len(components) != len(o.shape()):
            error("The number of components '%s' must be equal to the number of shapes '%s' for TensorConstant." % (str(components), str(o.shape())))

        # Let the UFL element handle the component map.
        component = o.element()._sub_element_mapping[components]

        # Handle restriction (offset by value shape).
        if self.restriction == "-":
            component += product(o.shape())

        # Let child class handle the constant coefficient
        return self.create_constant_coefficient(o.count(), component)

    # -------------------------------------------------------------------------
    # Indexed (indexed.py).
    # -------------------------------------------------------------------------
    def indexed(self, o):
        #print("\n\nVisiting Indexed:" + o.__repr__())

        # Get indexed expression and index, map index to current value and update components
        indexed_expr, index = o.operands()
        self._components.push(self.visit(index))

        #print "Indexed: indexed_expr: ", indexed_expr
        #print "Indexed: index: ", repr(index)
        #print "Indexed: comps: ", self._components

        # Visit expression subtrees and generate code.
        code = self.visit(indexed_expr)
        #print "Indexed: code: ", code

        # Remove component again
        self._components.pop()

        return code

    # -------------------------------------------------------------------------
    # MultiIndex (indexing.py).
    # -------------------------------------------------------------------------
    def multi_index(self, o):
        #print("\n\nVisiting MultiIndex:" + o.__repr__())

        # Loop all indices in MultiIndex and get current values
        subcomp = []
        for i in o:
            if isinstance(i, FixedIndex):
                subcomp.append(i._value)
            elif isinstance(i, Index):
                subcomp.append(self._index2value[i])
        #print "MultiIndex: subcomp: ", tuple(subcomp)
        return tuple(subcomp)

    # -------------------------------------------------------------------------
    # IndexSum (indexsum.py).
    # -------------------------------------------------------------------------
    def index_sum(self, o):
        #print("\n\nVisiting IndexSum: " + str(tree_format(o)))

        # Get expression and index that we're summing over
        summand, multiindex = o.operands()
        index, = multiindex

        #print "\nIndexSum: summand: ", summand
        #print "IndexSum: multiind: ", repr(multiindex)
        #print "IndexSum: index: ", repr(index)
        #print "IndexSum: o.dim: ", o.dimension()

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
        #print("\n\nVisiting Sqrt: " + o.__repr__() + "with operands: " + "\n".join(map(str,operands)))
        # Call common math function.
        return self._math_function(operands, self.format["sqrt"])

    def exp(self, o, *operands):
        #print("\n\nVisiting Exp: " + o.__repr__() + "with operands: " + "\n".join(map(str,operands)))
        # Call common math function.
        return self._math_function(operands, self.format["exp"])

    def ln(self, o, *operands):
        #print("\n\nVisiting Ln: " + o.__repr__() + "with operands: " + "\n".join(map(str,operands)))
        # Call common math function.
        return self._math_function(operands, self.format["ln"])

    def cos(self, o, *operands):
        #print("\n\nVisiting Cos: " + o.__repr__() + "with operands: " + "\n".join(map(str,operands)))
        # Call common math function.
        return self._math_function(operands, self.format["cos"])

    def sin(self, o, *operands):
        #print("\n\nVisiting Sin: " + o.__repr__() + "with operands: " + "\n".join(map(str,operands)))
        # Call common math function.
        return self._math_function(operands, self.format["sin"])

    # -------------------------------------------------------------------------
    # PositiveRestricted and NegativeRestricted (restriction.py).
    # -------------------------------------------------------------------------
    def positive_restricted(self, o):
        #print("\n\nVisiting PositiveRestricted: " + o.__repr__())

        # Just get the first operand, there should only be one.
        restricted_expr = o.operands()
        if len(restricted_expr) != 1:
            error("Only expected one operand for restriction: " + str(restricted_expr))
        if not self.restriction is None:
            error("Expression is restricted twice: " + str(restricted_expr))

        #print "PositiveRestricted expr: ", restricted_expr
 
        # Visit operand and generate restricted code.
        self.restriction = "+"
        code = self.visit(restricted_expr[0])
        #print "PositiveRestricted code: ", code

        # Reset restriction
        self.restriction = None

        return code

    def negative_restricted(self, o):
        #print("\n\nVisiting NegativeRestricted: " + o.__repr__())

        # Just get the first operand, there should only be one.
        restricted_expr = o.operands()
        if len(restricted_expr) != 1:
            error("Only expected one operand for restriction: " + str(restricted_expr))
 
        if not self.restriction is None:
            error("Expression is restricted twice: " + str(restricted_expr))

        # Visit operand and generate restricted code.
        self.restriction = "-"
        code = self.visit(restricted_expr[0])

        # Reset restriction
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

        #print "\nComponentTensor: component_expr: ", component_expr
        #print "ComponentTensor: indices: ", repr(indices)
        #print "ComponentTensor: components: ", components

        if len(components) != len(indices):
            error("The number of known components must be equal to the number of components of the ComponentTensor for this to work.")

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
        #print("\n\nVisiting ListTensor: " + o.__repr__())

        # Get the component
        component = self.component()
        #print "ListTensor: component: ", component

        # Extract first and the rest of the components
        c0, c1 = component[0], component[1:]
        #print "ListTensor: c0 ", c0
        #print "ListTensor: c1 ", c1

        # Get first operand
        op = o.operands()[c0]
        #print "ListTensor: op ", op

        # Evaluate subtensor with this subcomponent
        self._components.push(c1)
        code = self.visit(op)
        self._components.pop()

        return code

    # -------------------------------------------------------------------------
    # Variable (variable.py).
    # -------------------------------------------------------------------------
    def variable(self, o):
        #print("\n\nVisiting Variable: " + o.__repr__())
        # Just get the expression associated with the variable
        return self.visit(o.expression())

    # -------------------------------------------------------------------------
    # Functions that generates code from integrand
    # -------------------------------------------------------------------------
    def generate_code(self, integrand, Indent, interior):
        "Generate code from integrand."

        # Prefetch formats to speed up code generation.
        format_comment      = self.format["comment"]
        format_float_decl   = self.format["float declaration"]
        format_F            = self.format["function value"]
        format_float        = self.format["floating point"]
        format_add_equal    = self.format["add equal"]
        format_nzc          = self.format["nonzero columns"](0).split("0")[0]
        format_r            = self.format["free secondary indices"][0]
        format_mult         = self.format["multiply"]
        format_scale_factor = self.format["scale factor"]
        format_add          = self.format["add"]
        format_tensor       = self.format["element tensor quad"]
        format_array_access = self.format["array access"]
        format_Gip          = self.format["geometry tensor"] + self.format["integration points"]

        # Initialise return values.
        code = []
        num_ops = 0

        # Only propagate restrictions if we have an interior integral.
        if interior:
            integrand = propagate_restrictions(integrand)

        # Profiling
#        name = "test.prof"
#        prof = hotshot.Profile(name)
#        prof.runcall(self.visit, integrand)
#        prof.close()
#        stats = hotshot.stats.load(name)
##        stats.strip_dirs()
#        stats.sort_stats("time").print_stats(50)
#        raise RuntimeError

        # Generate loop code by transforming integrand.
        info("Transforming UFL integrand...")
        t = time.time()
        loop_code = self.visit(integrand)
        info("done, time = %f" % (time.time() - t))

        # Generate code.
        info("Generate code...")
        t = time.time()

        # TODO: Verify that test and trial functions will ALWAYS be rearranged to 0 and 1.
        indices = {-2: self.format["first free index"], -1: self.format["second free index"],
                    0: self.format["first free index"],  1: self.format["second free index"]}

        # Create the function declarations, we know that the code generator numbers
        # functions from 0 to n.
        if self.function_count:
            code += ["", format_comment("Function declarations")]
        for function_number in range(self.function_count):
            code.append((format_float_decl + format_F + str(function_number), format_float(0)))

        # Create code for computing function values, sort after loop ranges first.
        functions = self.functions
        function_list = {}
        for key, val in functions.items():
            if val[1] in function_list:
                function_list[val[1]].append(key)
            else:
                function_list[val[1]] = [key]

        # Loop ranges and get list of functions.
        for loop_range, list_of_functions in function_list.items():
            function_expr = {}
            function_numbers = []
            func_ops = 0
            # Loop functions.
            for function in list_of_functions:
                # Get name and number.
                name = str(functions[function][0])
                number = int(name.strip(format_F))
                # TODO: This check can be removed for speed later.
                if number in function_numbers:
                    error("This is definitely not supposed to happen!")
                function_numbers.append(number)
                # Get number of operations to compute entry and add to function operations count.
                f_ops = self._count_operations(function) + 1
                func_ops += f_ops
                entry = format_add_equal(name, function)
                function_expr[number] = entry

                # Extract non-zero column number if needed.
                if format_nzc in entry:
                    self.used_nzcs.add(int(entry.split(format_nzc)[1].split("[")[0]))

            # Multiply number of operations by the range of the loop index and add
            # number of operations to compute function values to total count.
            func_ops *= loop_range
            func_ops_comment = ["", format_comment("Total number of operations to compute function values = %d" % func_ops)]
            num_ops += func_ops

            # Sort the functions according to name and create loop to compute the function values.
            function_numbers.sort()
            lines = []
            for number in function_numbers:
                lines.append(function_expr[number])
            code += func_ops_comment + generate_loop(lines, [(format_r, 0, loop_range)], Indent, self.format)

        # Create weight.
        weight = self._weight()

        # Generate entries, multiply by weights and sort after primary loops.
        loops = {}
        for key, val in loop_code.items():

            # If value was zero continue.
            if val is None:
                continue

            # Create value, get number of operations and an indicator of zero valued value
            value, zero = self._create_value(val, weight, format_scale_factor)

            if zero:
                continue

            # Add points and scale factor to used weights and transformations
            self.used_weights.add(self.points)
            self.trans_set.add(format_scale_factor)

            # Compute number of operations to compute entry
            # (add 1 because of += in assignment).
            entry_ops = self._count_operations(value) + 1

            # Create comment for number of operations
            entry_ops_comment = format_comment("Number of operations to compute entry: %d" % entry_ops)

            # Create appropriate entries.
            # FIXME: We only support rank 0, 1 and 2.
            entry = ""
            loop = ()
            if len(key) == 0:
                entry = "0"

            elif len(key) == 1:
                key = key[0]
                # Checking if the basis was a test function.
                # TODO: Make sure test function indices are always rearranged to 0.
                if key[0] != -2 and key[0] != 0:
                    error("Linear forms must be defined using test functions only: " + str(key))

                index_j, entry, range_j, space_dim_j = key
                loop = ((indices[index_j], 0, range_j),)
                if range_j == 1 and self.optimise_options["ignore ones"]:
                    loop = ()
                # Multiply number of operations to compute entries by range of loop.
                entry_ops *= range_j

                # Extract non-zero column number if needed.
                if format_nzc in entry:
                    self.used_nzcs.add(int(entry.split(format_nzc)[1].split("[")[0]))

            elif len(key) == 2:
                # Extract test and trial loops in correct order and check if for is legal.
                key0, key1 = (0, 0)
                for k in key:
                    if not k[0] in indices:
                        error("Bilinear forms must be defined using test and trial functions (index -2, -1, 0, 1): " + str(k))
                    if k[0] == -2 or k[0] == 0:
                        key0 = k
                    else:
                        key1 = k
                index_j, entry_j, range_j, space_dim_j = key0
                index_k, entry_k, range_k, space_dim_k = key1

                loop = []
                if not (range_j == 1 and self.optimise_options["ignore ones"]):
                    loop.append((indices[index_j], 0, range_j))
                if not (range_k == 1 and self.optimise_options["ignore ones"]):
                    loop.append((indices[index_k], 0, range_k))

                entry = format_add([format_mult([entry_j, str(space_dim_k)]), entry_k])
                loop = tuple(loop)

                # Multiply number of operations to compute entries by range of loops.
                entry_ops *= range_j*range_k

                # Extract non-zero column number if needed.
                if format_nzc in entry_j:
                    self.used_nzcs.add(int(entry_j.split(format_nzc)[1].split("[")[0]))
                if format_nzc in entry_k:
                    self.used_nzcs.add(int(entry_k.split(format_nzc)[1].split("[")[0]))
            else:
                error("Only rank 0, 1 and 2 tensors are currently supported: " + str(key))

            # Generate the code line for the entry.
            # Try to evaluate entry ("3*6 + 2" --> "20").
            try:
                entry = str(eval(entry))
            except:
                pass

            entry_code = format_add_equal( format_tensor + format_array_access(entry), value)

            if loop not in loops:
                loops[loop] = [entry_ops, [entry_ops_comment, entry_code]]
            else:
                loops[loop][0] += entry_ops
                loops[loop][1] += [entry_ops_comment, entry_code]

        # Generate code for ip constant declarations.
        ip_const_ops, ip_const_code = generate_aux_constants(self.ip_consts, format_Gip,\
                                        self.format["const float declaration"], True)
        num_ops += ip_const_ops
        if ip_const_code:
            code += ["", format_comment("Number of operations to compute ip constants: %d" %ip_const_ops)]
            code += ip_const_code

        # Write all the loops of basis functions.
        for loop, ops_lines in loops.items():
            ops, lines = ops_lines

            # Add number of operations for current loop to total count.
            num_ops += ops
            code += ["", format_comment("Number of operations for primary indices: %d" % ops)]
            code += generate_loop(lines, loop, Indent, self.format)

        info("             done, time = %f" % (time.time() - t))

        # Reset ip constant declarations
        self.ip_consts = {}

        return code, num_ops

    def _count_operations(self, expression):
        error("This function should be implemented by the child class.")

    def _weight(self):
        error("This function should be implemented by the child class.")

    def _create_value(self, val, weight, scale_factor):
        error("This function should be implemented by the child class.")

