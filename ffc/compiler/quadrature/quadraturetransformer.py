"QuadratureTransformer for quadrature code generation to translate UFL expressions."

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2009-02-09 -- 2009-08-08"
__copyright__ = "Copyright (C) 2009 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Python modules.
from numpy import shape

# UFL common.
from ufl.common import product

# UFL Classes.
from ufl.classes import MultiIndex, FixedIndex, IntValue, FloatValue, Function

# UFL Algorithms.
from ufl.algorithms.transformations import Transformer
from ufl.algorithms import purge_list_tensors, expand_indices, propagate_restrictions, strip_variables
from ufl.algorithms.printing import tree_format

# FFC common modules.
from ffc.common.log import info, debug, error

# FFC compiler modules.
from ffc.compiler.tensor.multiindex import MultiIndex as FFCMultiIndex

# FFC fem modules.
from ffc.fem.createelement import create_element
from ffc.fem.finiteelement import AFFINE, CONTRAVARIANT_PIOLA, COVARIANT_PIOLA

# Utility and optimisation functions for quadraturegenerator.
from quadraturegenerator_utils import generate_loop, generate_psi_name, create_psi_tables, create_permutations
from reduce_operations import operation_count

import time

class QuadratureTransformer(Transformer):
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
    # Supported classes.
    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------
    # AlgebraOperators (algebra.py).
    # -------------------------------------------------------------------------
    def sum(self, o, *operands):
        debug("Visiting Sum: " + o.__repr__() + "\noperands: " + "\n".join(map(str, operands)))

        # Prefetch formats to speed up code generation.
        format_group  = self.format["grouping"]
        format_add    = self.format["add"]
        code = {}

        # Loop operands that has to be summend.
        for op in operands:
            # If entries does already exist we can add the code, otherwise just
            # dump them in the element tensor.
            for key, val in op.items():
                if key in code:
                    code[key].append(val)
                else:
                    code[key] = [val]

        # Add sums and group if necessary.
        for key, val in code.items():
            if len(val) > 1:
                value = [v for v in val if v]
                code[key] = format_group(format_add(value))
            else:
                code[key] = val[0]

        return code

    def product(self, o, *operands):
        debug("\n\nVisiting Product: " + o.__repr__() + "with operands: " + "\n".join(map(str,operands)))

        # Prefetch formats to speed up code generation.
        format_mult = self.format["multiply"]
        permute = []
        not_permute = []

        # Sort operands in objects that needs permutation and objects that does not.
        for op in operands:
            if len(op) > 1 or (op and op.keys()[0] != ()):
                permute.append(op)
            elif op:
                not_permute.append(op[()])

        # Create permutations.
        permutations = create_permutations(permute)

        debug("\npermute: " + str(permute))
        debug("\nnot_permute: " + str(not_permute))
        debug("\npermutations: " + str(permutations))

        # Create code.
        code ={}
        if permutations:
            for key, val in permutations.items():
                # Sort key in order to create a unique key.
                l = list(key)
                l.sort()
                # TODO: I think this check can be removed for speed since we
                # just have a list of objects we should never get any conflicts here.
                value = [v for v in val + not_permute if v]
                if tuple(l) in code:
                    error("This key should not be in the code.")
                code[tuple(l)] = format_mult(value)
        else:
            value = [v for v in not_permute if v]
            code[()] = format_mult(value)

        return code

    def division(self, o, *operands):
        debug("\n\nVisiting Division: " + o.__repr__() + "with operands: " + "\n".join(map(str,operands)))

        # Prefetch formats to speed up code generation.
        format_div      = self.format["division"]
        format_grouping = self.format["grouping"]

        if len(operands) != 2:
            error("Expected exactly two operands (numerator and denominator): " + operands.__repr__())

        # Get the code from the operands.
        numerator_code, denominator_code = operands
        debug("\nnumerator: " + str(numerator_code))
        debug("\ndenominator: " + str(denominator_code))

        # TODO: Are these safety checks needed?
        if not () in denominator_code and len(denominator_code) != 1:
            error("Only support function type denominator: " + str(denominator_code))

        # Get denominator and create new values for the numerator.
        denominator = denominator_code.pop(())
        for key, val in numerator_code.items():
            numerator_code[key] = val + format_div + format_grouping(denominator)

        return numerator_code

    def power(self, o):
        debug("\n\nVisiting Power: " + o.__repr__())

        # Get base and exponent.
        base, expo = o.operands()
        debug("\nbase: " + str(base))
        debug("\nexponent: " + str(expo))

        # Visit base to get base code.
        base_code = self.visit(base)
        debug("base_code: " + str(base_code))

        # TODO: Are these safety checks needed?
        if not () in base_code and len(base_code) != 1:
            error("Only support function type base: " + str(base_code))

        # Get the base code.
        val = base_code.pop(())

        # Handle different exponents
        if isinstance(expo, IntValue):
            return {(): self.format["power"](val, expo.value())}
        elif isinstance(expo, FloatValue):
            return {(): self.format["std power"](val, self.format["floating point"](expo.value()))}
        elif isinstance(expo, Function):
            exp = self.visit(expo)
            return {(): self.format["std power"](val, exp.pop(()))}
        else:
            error("power does not support this exponent: " + repr(expo))

    def abs(self, o, *operands):
        debug("\n\nVisiting Abs: " + o.__repr__() + "with operands: " + "\n".join(map(str,operands)))

        # Prefetch formats to speed up code generation.
        format_abs = self.format["absolute value"]

        # TODO: Are these safety checks needed?
        if len(operands) != 1 and not () in operands[0] and len(operands[0]) != 1:
            error("Abs expects one operand of function type: " + str(operands))

        # Take absolute value of operand.
        operand = operands[0]
        for key, val in operand.items():
            operand[key] = format_abs(val)

        return operand

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
    # Constant values (constantvalue.py).
    # -------------------------------------------------------------------------
    def scalar_value(self, o, *operands):
        "ScalarValue covers IntValue and FloatValue"
        debug("\n\nVisiting FloatValue:" + o.__repr__())

        # FIXME: Might be needed because it can be IndexAnnotated?
        if operands:
            error("Did not expect any operands for ScalarValue: " + str((o, operands)))

        # TODO: Handle value < 0 better such that we don't have + -2 in the code.
        return {():self.format["floating point"](o.value())}

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
    # Function and Constants (function.py).
    # -------------------------------------------------------------------------
    def constant(self, o, *operands):
        debug("\n\nVisiting Constant: " + o.__repr__())

        # Safety checks.
        if operands:
            error("Didn't expect any operands for Constant: " + str(operands))
        if len(self._components) > 0:
            error("Constant does not expect component indices: " + str(self._components))
        if o.shape() != ():
            error("Constant should not have a value shape: " + str(o.shape()))

        component = 0
        # Handle restriction.
        if self.restriction == "-":
            component += 1

        coefficient = self.format["coeff"] + self.format["matrix access"](str(o.count()), component)
        debug("Constant coefficient: " + coefficient)
        return {():coefficient}

    def vector_constant(self, o, *operands):
        debug("\n\nVisiting VectorConstant: " + o.__repr__())
        # Safety checks.
        if operands:
            error("Didn't expect any operands for VectorConstant: " + str(operands))
        if len(self._components) != 1 or not isinstance(self._components[0], FixedIndex):
            error("VectorConstant expects 1 Fixed component index: " + str(self._components))

        # We get one component.
        component = int(self._components[0])

        # Handle restriction.
        if self.restriction == "-":
            component += o.shape()[0]

        coefficient = self.format["coeff"] + self.format["matrix access"](str(o.count()), component)
        debug("VectorConstant coefficient: " + coefficient)
        return {():coefficient}

    def tensor_constant(self, o, *operands):
        debug("\n\nVisiting TensorConstant: " + o.__repr__())

        # Safety checks.
        if operands:
            error("Didn't expect any operands for TensorConstant: " + str(operands))
        if not all(isinstance(c, FixedIndex) for c in self._components):
            error("TensorConstant expects FixedIndex as components: " + str(self._components))

        # Compute the global component.
        component = tuple([int(c) for c in self._components])
        component = o.element()._sub_element_mapping[component]

        # Handle restriction (offset by value shape).
        if self.restriction == "-":
            component += product(o.shape())

        coefficient = self.format["coeff"] + self.format["matrix access"](str(o.count()), component)
        debug("TensorConstant coefficient: " + coefficient)
        return {():coefficient}

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
    # FacetNormal (geometry.py).
    # -------------------------------------------------------------------------
    def facet_normal(self, o,  *operands):
        debug("Visiting FacetNormal:")
        # Safety checks.
        if operands:
            error("Didn't expect any operands for FacetNormal: " + str(operands))
        if len(self._components) != 1 or not isinstance(self._components[0], FixedIndex):
            error("FacetNormal expects 1 Fixed component index: " + str(self._components))

        # We get one component.
        component = int(self._components[0])
        normal_component = self.format["normal component"] + str(component)
        self.trans_set.add(normal_component)
        debug("Facet Normal Component: " + normal_component)
        return {():normal_component}

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
    def _math_function(self, operands, format_function):
        # TODO: Are these safety checks needed?
        if len(operands) != 1 and not () in operands[0] and len(operands[0]) != 1:
            error("MathFunctions expect one operand of function type: " + str(operands))
        # Use format function on value of operand.
        operand = operands[0]
        for key, val in operand.items():
            operand[key] = format_function(val)
        debug("operand: " + str(operand))
        return operand

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

    # -------------------------------------------------------------------------
    # Helper functions for BasisFunction and Function).
    # -------------------------------------------------------------------------
    def create_basis_function(self, ufl_basis_function, component, derivatives):
        "Create code for basis functions, and update relevant tables of used basis."

        # Prefetch formats to speed up code generation.
        format_group         = self.format["grouping"]
        format_add           = self.format["add"]
        format_mult          = self.format["multiply"]
        format_transform     = self.format["transform"]
        format_detJ          = self.format["determinant"]
        format_inv           = self.format["inverse"]

        # Get local component (in case we have mixed elements).
        local_comp, local_elem = ufl_basis_function.element().extract_component(tuple(component))

        # Check that we don't take derivatives of QuadratureElements.
        if derivatives and local_elem.family() == "Quadrature":
            error("Derivatives of Quadrature elements are not supported: " + str(ufl_basis_function))

        # Handle tensor elements.
        if len(local_comp) > 1:
            local_comp = local_elem._sub_element_mapping[local_comp]
        elif local_comp:
            local_comp = local_comp[0]
        else:
            local_comp = 0

        local_offset = 0
        if len(component) > 1:
            component = ufl_basis_function.element()._sub_element_mapping[tuple(component)]
        elif component:
            component = component.pop()

        # Compute the local offset, needed for non-affine mappings because the
        # elements are labeled with the global component number.
        if component:
            local_offset = component - local_comp

        # Create FFC element.
        ffc_element = create_element(ufl_basis_function.element())

        code = {}
        # Set geo_dim.
        # TODO: All terms REALLY have to be defined on cell with the same
        # geometrical dimension so only do this once and exclude the check?
        geo_dim = ufl_basis_function.element().cell().geometric_dimension()
        if self.geo_dim:
            if geo_dim != self.geo_dim:
                error("All terms must be defined on cells with the same geometrical dimension.")
        else:
            self.geo_dim = geo_dim

        # Generate FFC multi index for derivatives.
        multiindices = FFCMultiIndex([range(geo_dim)]*len(derivatives)).indices

        # Loop derivatives and get multi indices.
        for multi in multiindices:
            deriv = [multi.count(i) for i in range(geo_dim)]
            if not any(deriv):
                deriv = []
            transformation = ffc_element.component_element(component)[0].mapping()
            if transformation == AFFINE:
                # Call function to create mapping and basis name.
                mapping, basis = self.__create_mapping_basis(component, deriv, ufl_basis_function, ffc_element)
                if basis == None:
                    continue
                # Add transformation if needed.
                transforms = []
                for i, direction in enumerate(derivatives):
                    ref = multi[i]
                    t = format_transform("JINV", ref, direction, self.restriction)
                    self.trans_set.add(t)
                    transforms.append(t)
                # Only multiply by basis if it is present.
                if basis:
                    prods = transforms + [basis]
                else:
                    prods = transforms

                if mapping in code:
                    code[mapping].append(format_mult(prods))
                else:
                    code[mapping] = [format_mult(prods)]
            # Handle non-affine mappings.
            else:
                for c in range(geo_dim):
                    # Call function to create mapping and basis name.
                    mapping, basis = self.__create_mapping_basis(c + local_offset, deriv, ufl_basis_function, ffc_element)
                    if basis == None:
                        continue

                    # Multiply basis by appropriate transform.
                    if transformation == COVARIANT_PIOLA:
                        dxdX = format_transform("JINV", c, local_comp, self.restriction)
                        self.trans_set.add(dxdX)
                        basis = format_mult([dxdX, basis])
                    elif transformation == CONTRAVARIANT_PIOLA:
                        self.trans_set.add(format_detJ(self.restriction))
                        detJ = format_inv(format_detJ(self.restriction))
                        dXdx = format_transform("J", c, local_comp, self.restriction)
                        self.trans_set.add(dXdx)
                        basis = format_mult([detJ, dXdx, basis])
                    else:
                        error("Transformation is not supported: " + str(transformation))

                    # Add transformation if needed.
                    transforms = []
                    for i, direction in enumerate(derivatives):
                        ref = multi[i]
                        t = format_transform("JINV", ref, direction, self.restriction)
                        self.trans_set.add(t)
                        transforms.append(t)

                    # Only multiply by basis if it is present.
                    if basis:
                        prods = transforms + [basis]
                    else:
                        prods = transforms

                    if mapping in code:
                        code[mapping].append(format_mult(prods))
                    else:
                        code[mapping] = [format_mult(prods)]

        # Add sums and group if necessary.
        for key, val in code.items():
            if len(val) > 1:
                code[key] = format_group(format_add(val))
            else:
                code[key] = val[0]
        return code

    def __create_mapping_basis(self, component, deriv, ufl_basis_function, ffc_element):
        "Create basis name and mapping from given basis_info."

        # Get string for integration points.
        format_ip = self.format["integration points"]

        # Only support test and trial functions.
        # TODO: Verify that test and trial functions will ALWAYS be rearranged to 0 and 1.
        indices = {-2: self.format["first free index"],
                   -1: self.format["second free index"],
                    0: self.format["first free index"],
                    1: self.format["second free index"]}

        # Check that we have a basis function.
        if not ufl_basis_function.count() in indices:
            error("Currently, BasisFunction index must be either -2, -1, 0 or 1: " + str(ufl_basis_function))

        # Handle restriction through facet.
        facet = {"+": self.facet0, "-": self.facet1, None: self.facet0}[self.restriction]

        # Get element counter and loop index.
        element_counter = self.element_map[self.points][ufl_basis_function.element()]
        loop_index = indices[ufl_basis_function.count()]

        # Create basis access, we never need to map the entry in the basis table
        # since we will either loop the entire space dimension or the non-zeros.
        if self.points == 1:
            format_ip = "0"
        basis_access = self.format["matrix access"](format_ip, loop_index)

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

        if zeros and self.optimise_options["ignore zero tables"]:
            return (None, None)

        # If the loop index range is one we can look up the first component
        # in the psi array. If we only have ones we don't need the basis.
        basis = ""
        if self.optimise_options["ignore ones"] and loop_index_range == 1 and ones:
            loop_index = "0"
        else:
            # Add basis name to set of used tables and add matrix access.
            self.used_psi_tables.add(name)
            basis = name + basis_access

        # Create the correct mapping of the basis function into the local element tensor.
        basis_map = loop_index
        if non_zeros and basis_map == "0":
            basis_map = str(non_zeros[1][0])
        elif non_zeros:
            basis_map = self.format["nonzero columns"](non_zeros[0]) +\
                        self.format["array access"](basis_map)
        if offset:
            basis_map = self.format["grouping"](self.format["add"]([basis_map, offset]))

        # Try to evaluate basis map ("3 + 2" --> "5").
        try:
            basis_map = str(eval(basis_map))
        except:
            pass

        # Create mapping (index, map, loop_range, space_dim).
        # Example dx and ds: (0, j, 3, 3)
        # Example dS: (0, (j + 3), 3, 6), 6=2*space_dim
        # Example dS optimised: (0, (nz2[j] + 3), 2, 6), 6=2*space_dim
        mapping = ((ufl_basis_function.count(), basis_map, loop_index_range, space_dim),)

        return (mapping, basis)

    def create_function(self, ufl_function, component, derivatives):
        "Create code for basis functions, and update relevant tables of used basis."

        # Prefetch formats to speed up code generation.
        format_mult          = self.format["multiply"]
        format_transform     = self.format["transform"]
        format_detJ          = self.format["determinant"]
        format_inv           = self.format["inverse"]

        # Get local component (in case we have mixed elements).
        local_comp, local_elem = ufl_function.element().extract_component(tuple(component))

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

        local_offset = 0
        if len(component) > 1:
            component = ufl_function.element()._sub_element_mapping[tuple(component)]
        elif component:
            component = component.pop()

        # Compute the local offset (needed for non-affine mappings).
        if component:
            local_offset = component - local_comp

        # Create FFC element.
        ffc_element = create_element(ufl_function.element())
        code = []

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
        for multi in multiindices:
            deriv = [multi.count(i) for i in range(geo_dim)]
            if not any(deriv):
                deriv = []
            transformation = ffc_element.component_element(component)[0].mapping()
            if transformation == AFFINE:
                # Call other function to create function name.
                function_name = self.__create_function_name(component, deriv, quad_element, ufl_function, ffc_element)
                if not function_name:
                    continue

                # Add transformation if needed.
                transforms = []
                for i, direction in enumerate(derivatives):
                    ref = multi[i]
                    t = format_transform("JINV", ref, direction, self.restriction)
                    self.trans_set.add(t)
                    transforms.append(t)

                # Multiply function value by the transformations and add to code.
                code.append(format_mult(transforms + [function_name]))

            # Handle non-affine mappings.
            else:
                for c in range(geo_dim):
                    function_name = self.__create_function_name(c + local_offset, deriv, quad_element, ufl_function, ffc_element)

                    # Multiply basis by appropriate transform.
                    if transformation == COVARIANT_PIOLA:
                        dxdX = format_transform("JINV", c, local_comp, self.restriction)
                        self.trans_set.add(dxdX)
                        function_name = format_mult([dxdX, function_name])
                    elif transformation == CONTRAVARIANT_PIOLA:
                        self.trans_set.add(format_detJ(self.restriction))
                        detJ = format_inv(format_detJ(self.restriction))
                        dXdx = format_transform("J", c, local_comp, self.restriction)
                        self.trans_set.add(dXdx)
                        function_name = format_mult([detJ, dXdx, function_name])
                    else:
                        error("Transformation is not supported: ", str(transformation))

                    # Add transformation if needed.
                    transforms = []
                    for i, direction in enumerate(derivatives):
                        ref = multi[i]
                        t = format_transform("JINV", ref, direction, self.restriction)
                        self.trans_set.add(t)
                        transforms.append(t)

                    # Multiply function value by the transformations and add to code.
                    code.append(format_mult(transforms + [function_name]))

        if not code:
            return "0"
        elif len(code) > 1:
            code = self.format["grouping"](self.format["add"](code))
        else:
            code = code[0]

        return code

    def __create_function_name(self, component, deriv, quad_element, ufl_function, ffc_element):

        # Get string for integration points.
        format_ip = self.format["integration points"]

        # Pick first free index of secondary type
        # (could use primary indices, but it's better to avoid confusion).
        loop_index = self.format["free secondary indices"][0]

        # Create basis access, we never need to map the entry in the basis
        # table since we will either loop the entire space dimension or the
        # non-zeros.
        if self.points == 1:
            format_ip = "0"
        basis_access = self.format["matrix access"](format_ip, loop_index)

        # Handle restriction through facet.
        facet = {"+": self.facet0, "-": self.facet1, None: self.facet0}[self.restriction]

        # Get the element counter.
        element_counter = self.element_map[self.points][ufl_function.element()]

        # Offset by element space dimension in case of negative restriction.
        offset = {"+": "", "-": str(ffc_element.space_dimension()), None: ""}[self.restriction]

        # Create basis name and map to correct basis and get info.
        basis_name = generate_psi_name(element_counter, facet, component, deriv)
        basis_name, non_zeros, zeros, ones = self.name_map[basis_name]

        # If all basis are zero we just return "0".
        # TODO: Handle this more elegantly such that all terms involving this
        # zero factor is removed.
        if zeros and self.optimise_options["ignore zero tables"]:
            return None

        # Get the index range of the loop index.
        loop_index_range = shape(self.unique_tables[basis_name])[1]

        # Set default coefficient access.
        coefficient_access = loop_index

        # If the loop index range is one we can look up the first component
        # in the coefficient array. If we only have ones we don't need the basis.
        if self.optimise_options["ignore ones"] and loop_index_range == 1 and ones:
            coefficient_access = "0"
            basis_name = ""
        elif not quad_element:
            # Add basis name to set of used tables and add matrix access.
            self.used_psi_tables.add(basis_name)
            basis_name += basis_access

        # If we have a quadrature element we can use the ip number to look
        # up the value directly. Need to add offset in case of components.
        if quad_element:
            quad_offset = 0
            if component:
                for i in range(component):
                    quad_offset += ffc_element.sub_element(i).space_dimension()
            if quad_offset:
                coefficient_access = self.format["add"]([format_ip, str(quad_offset)])
            else:
                coefficient_access = format_ip

        # If we have non zero column mapping but only one value just pick it.
        if non_zeros and coefficient_access == "0":
            coefficient_access = str(non_zeros[1][0])
        elif non_zeros and not quad_element:
            coefficient_access = self.format["nonzero columns"](non_zeros[0]) +\
                                 self.format["array access"](coefficient_access)
        if offset:
            coefficient_access = self.format["add"]([coefficient_access, offset])

        # Try to evaluate coefficient access ("3 + 2" --> "5").
        try:
            coefficient_access = str(eval(coefficient_access))
        except:
            pass

        coefficient = self.format["coeff"] +\
                      self.format["matrix access"](str(ufl_function.count()), coefficient_access)
        function_expr = coefficient
        if basis_name:
            function_expr = self.format["multiply"]([basis_name, coefficient])

        # If we have a quadrature element (or if basis was deleted) we don't need the basis.
        if quad_element or not basis_name:
            function_name = coefficient
        else:
            # Check if the expression to compute the function value is already in
            # the dictionary of used function. If not, generate a new name and add.
            function_name = self.format["function value"] + str(self.function_count)
            if not function_expr in self.functions:
                self.functions[function_expr] = (function_name, loop_index_range)
                # Increase count.
                self.function_count += 1
            else:
                function_name, index_r = self.functions[function_expr]
                # Check just to make sure.
                if not index_r == loop_index_range:
                    error("Index ranges does not match.")
        return function_name

def generate_code(integrand, transformer, Indent, format, interior):
    """Generate code from a UFL integral type. It generates all the code that
    goes inside the quadrature loop."""

    # Prefetch formats to speed up code generation.
    format_comment      = format["comment"]
    format_float_decl   = format["float declaration"]
    format_F            = format["function value"]
    format_float        = format["floating point"]
    format_add_equal    = format["add equal"]
    format_nzc          = format["nonzero columns"](0).split("0")[0]
    format_r            = format["free secondary indices"][0]
    format_mult         = format["multiply"]
    format_scale_factor = format["scale factor"]
    format_add          = format["add"]
    format_tensor       = format["element tensor quad"]
    format_array_access = format["array access"]

    # Initialise return values.
    code = []
    num_ops = 0

    debug("\nQG, Using Transformer.")

    # Apply basic expansions.
    # TODO: Figure out if there is a 'correct' order of doing this
    # In form.form_data().form, which we should be using, coefficients have
    # been mapped and derivatives expandes. So it should be enough to just
    # expand_indices and purge_list_tensors.
#    t = time.time()
    new_integrand = strip_variables(integrand)
    new_integrand = expand_indices(new_integrand)
#    info("expand_indices, time = %f" % (time.time() - t))
#    t = time.time()
    new_integrand = purge_list_tensors(new_integrand)
#    info("purge_tensors, time  = %f" % (time.time() - t))
    # Only propagate restrictions if we have an interior integral.
    if interior:
        new_integrand = propagate_restrictions(new_integrand)
    debug("\nExpanded integrand\n" + str(tree_format(new_integrand)))

    # Let the Transformer create the loop code.
#    t = time.time()
    loop_code = transformer.visit(new_integrand)
#    info("gen. loop_code, time = %f" % (time.time() - t))

    # TODO: Verify that test and trial functions will ALWAYS be rearranged to 0 and 1.
    indices = {-2: format["first free index"], -1: format["second free index"],
                0: format["first free index"],  1: format["second free index"]}

    # Create the function declarations, we know that the code generator numbers
    # functions from 0 to n.
    if transformer.function_count:
        code += ["", format_comment("Function declarations")]
    for function_number in range(transformer.function_count):
        code.append((format_float_decl + format_F + str(function_number), format_float(0)))

    # Create code for computing function values, sort after loop ranges first.
    functions = transformer.functions
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
            name = functions[function][0]
            number = int(name.strip(format_F))
            # TODO: This check can be removed for speed later.
            if number in function_numbers:
                error("This is definitely not supposed to happen!")
            function_numbers.append(number)
            # Get number of operations to compute entry and add to function operations count.
            f_ops = operation_count(function, format) + 1
            func_ops += f_ops
            entry = format_add_equal(name, function)
            function_expr[number] = entry

            # Extract non-zero column number if needed.
            if format_nzc in entry:
                transformer.used_nzcs.add(int(entry.split(format_nzc)[1].split("[")[0]))

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
        code += func_ops_comment + generate_loop(lines, [(format_r, 0, loop_range)], Indent, format)

    # Create weight.
    weight = format["weight"](transformer.points)
    if transformer.points > 1:
        weight += format_array_access(format["integration points"])

    # Generate entries, multiply by weights and sort after primary loops.
    loops = {}
#    t = time.time()
    for key, val in loop_code.items():

        # If value was zero continue.
        if val == None:
            continue
        # Multiply by weight and determinant, add both to set of used weights and transforms.
        value = format_mult([val, weight, format_scale_factor])
        transformer.used_weights.add(transformer.points)
        transformer.trans_set.add(format_scale_factor)

        # Compute number of operations to compute entry and create comment
        # (add 1 because of += in assignment).
        entry_ops = operation_count(value, format) + 1
        entry_ops_comment = format_comment("Number of operations to compute entry = %d" % entry_ops)
        prim_ops = entry_ops

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
            if range_j == 1 and transformer.optimise_options["ignore ones"]:
                loop = ()
            # Multiply number of operations to compute entries by range of loop.
            prim_ops *= range_j

            # Extract non-zero column number if needed.
            if format_nzc in entry:
                transformer.used_nzcs.add(int(entry.split(format_nzc)[1].split("[")[0]))

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
            if not (range_j == 1 and transformer.optimise_options["ignore ones"]):
                loop.append((indices[index_j], 0, range_j))
            if not (range_k == 1 and transformer.optimise_options["ignore ones"]):
                loop.append((indices[index_k], 0, range_k))

            entry = format_add([format_mult([entry_j, str(space_dim_k)]), entry_k])
            loop = tuple(loop)

            # Multiply number of operations to compute entries by range of loops.
            prim_ops *= range_j*range_k

            # Extract non-zero column number if needed.
            if format_nzc in entry_j:
                transformer.used_nzcs.add(int(entry_j.split(format_nzc)[1].split("[")[0]))
            if format_nzc in entry_k:
                transformer.used_nzcs.add(int(entry_k.split(format_nzc)[1].split("[")[0]))
        else:
            error("Only rank 0, 1 and 2 tensors are currently supported: " + str(key))

        # Generate the code line for the entry.
        entry_code = format_add_equal( format_tensor + format_array_access(entry), value)

        if loop not in loops:
            loops[loop] = [prim_ops, [entry_ops_comment, entry_code]]
        else:
            loops[loop][0] += prim_ops
            loops[loop][1] += [entry_ops_comment, entry_code]

    # Write all the loops of basis functions.
    for loop, ops_lines in loops.items():
        ops, lines = ops_lines

        # Add number of operations for current loop to total count.
        num_ops += ops
        code += ["", format_comment("Number of operations for primary indices = %d" % ops)]
        code += generate_loop(lines, loop, Indent, format)
#    info("write code, time     = %f" % (time.time() - t))
    return (code, num_ops)

