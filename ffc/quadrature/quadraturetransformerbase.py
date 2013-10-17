"""QuadratureTransformerBase, a common class for quadrature
transformers to translate UFL expressions."""

# Copyright (C) 2009-2013 Kristian B. Oelgaard
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Martin Alnaes, 2013
# Modified by Garth N. Wells, 2013
#
# First added:  2009-10-13
# Last changed: 2013-02-20

# Python modules.
from itertools import izip
from numpy import shape, array

# UFL Classes.
from ufl.classes import FixedIndex, Index
from ufl.common import StackDict, Stack, product
from ufl.permutation import build_component_numbering

# UFL Algorithms.
from ufl.algorithms import Transformer

# FFC modules.
from ffc.log import ffc_assert, error, info
from ffc.fiatinterface import create_element, map_facet_points
from ffc.mixedelement import MixedElement
from ffc.cpp import format

# FFC tensor modules.
from ffc.tensor.multiindex import MultiIndex as FFCMultiIndex
from ffc.representationutils import transform_component

# Utility and optimisation functions for quadraturegenerator.
from quadratureutils import create_psi_tables
from symbolics import BASIS, IP, GEO, CONST

class QuadratureTransformerBase(Transformer):
    "Transform UFL representation to quadrature code."

    def __init__(self,
                 psi_tables,
                 quad_weights,
                 gdim,
                 tdim,
                 entitytype,
                 function_replace_map,
                 optimise_parameters):

        Transformer.__init__(self)

        # Save optimise_parameters, weights and fiat_elements_map.
        self.optimise_parameters = optimise_parameters

        # Map from original functions with possibly incomplete elements
        # to functions with properly completed elements
        self._function_replace_map = function_replace_map
        self._function_replace_values = set(function_replace_map.values()) # For assertions

        # Create containers and variables.
        self.used_psi_tables = set()
        self.psi_tables_map = {}
        self.used_weights = set()
        self.quad_weights = quad_weights
        self.used_nzcs = set()
        self.ip_consts = {}
        self.trans_set = set()
        self.function_data = {}
        self.tdim = tdim
        self.gdim = gdim
        self.entitytype = entitytype
        self.points = 0
        self.facet0 = None
        self.facet1 = None
        self.vertex = None
        self.restriction = None
        self.avg = None
        self.coordinate = None
        self.conditionals = {}
        self.additional_includes_set = set()
        self.__psi_tables = psi_tables # TODO: Unused? Remove?

        # Stacks.
        self._derivatives = []
        self._index2value = StackDict()
        self._components = Stack()

        self.element_map, self.name_map, self.unique_tables =\
            create_psi_tables(psi_tables, self.optimise_parameters["eliminate zeros"], self.entitytype)

        # Cache.
        self.argument_cache = {}
        self.function_cache = {}

    def update_cell(self):
        ffc_assert(self.entitytype == "cell", "Not expecting update_cell on a %s." % self.entitytype)
        self.facet0 = None
        self.facet1 = None
        self.vertex = None
        self.coordinate = None
        self.conditionals = {}

    def update_facets(self, facet0, facet1):
        ffc_assert(self.entitytype == "facet", "Not expecting update_facet on a %s." % self.entitytype)
        self.facet0 = facet0
        self.facet1 = facet1
        self.vertex = None
        self.coordinate = None
        self.conditionals = {}

    def update_vertex(self, vertex):
        ffc_assert(self.entitytype == "vertex", "Not expecting update_vertex on a %s." % self.entitytype)
        self.facet0 = None
        self.facet1 = None
        self.vertex = vertex
        self.coordinate = None
        self.conditionals = {}

    def update_points(self, points):
        self.points = points
        self.coordinate = None
        # Reset functions everytime we move to a new quadrature loop
        self.conditionals = {}
        self.function_data = {}

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
        error("This expression is not handled: " + repr(o))

    # Nothing in terminal.py is handled. Can only handle children of these clases.
    def terminal(self, o):
        print "\n\nVisiting basic Terminal:", repr(o), "with operands:"
        error("This terminal is not handled: " + repr(o))

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
        error("All derivatives apart from Grad should have been expanded!!")

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
        error("This type of Condition is not supported (yet).")

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

    def math_function(self, o):
        print "\n\nVisiting MathFunction:", repr(o)
        error("This MathFunction is not supported (yet).")

    def atan_2_function(self, o):
        print "\n\nVisiting Atan2Function:", repr(o)
        error("Atan2Function is not implemented (yet).")

    def bessel_function(self, o):
        print "\n\nVisiting BesselFunction:", repr(o)
        error("BesselFunction is not implemented (yet).")

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
    # FacetNormal, CellVolume, Circumradius (geometry.py).
    # -------------------------------------------------------------------------
    def facet_normal(self, o):
        print "\n\nVisiting FacetNormal: ", repr(o)
        error("This object should be implemented by the child class.")

    def cell_volume(self, o):
        print "\n\nVisiting CellVolume: ", repr(o)
        error("This object should be implemented by the child class.")

    def circumradius(self, o):
        print "\n\nVisiting Circumeradius: ", repr(o)
        error("This object should be implemented by the child class.")

    # -------------------------------------------------------------------------
    # Things that can be handled by the base class.
    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------
    # Argument (basisfunction.py).
    # -------------------------------------------------------------------------
    def argument(self, o):
        #print("\nVisiting Argument:" + repr(o))

        # Map o to object with proper element and numbering
        o = self._function_replace_map[o]

        # Create aux. info.
        components = self.component()
        derivatives = self.derivatives()

        # Check if basis is already in cache
        key = (o, components, derivatives, self.restriction, self.avg)
        basis = self.argument_cache.get(key, None)

        # FIXME: Why does using a code dict from cache make the expression manipulations blow (MemoryError) up later?
        if basis is None or self.optimise_parameters["optimisation"]:
            # Get auxiliary variables to generate basis
            (component, local_elem, local_comp, local_offset,
             ffc_element, transformation, multiindices) = self._get_auxiliary_variables(o, components, derivatives)

            # Create mapping and code for basis function and add to dict.
            basis = self.create_argument(o, derivatives, component, local_comp,
                                         local_offset, ffc_element,
                                         transformation, multiindices,
                                         self.tdim, self.gdim, self.avg)
            self.argument_cache[key] = basis

        return basis

    # -------------------------------------------------------------------------
    # Constant values (constantvalue.py).
    # -------------------------------------------------------------------------
    def identity(self, o):
        #print "\n\nVisiting Identity: ", repr(o)

        # Get components
        i, j = self.component()

        # Only return a value if i==j
        if i == j:
            return self._format_scalar_value(1.0)
        else:
            return self._format_scalar_value(None)

    def scalar_value(self, o):
        "ScalarValue covers IntValue and FloatValue"
        #print "\n\nVisiting ScalarValue: ", repr(o)
        return self._format_scalar_value(o.value())

    def zero(self, o):
        #print "\n\nVisiting Zero:", repr(o)
        return self._format_scalar_value(None)

    # -------------------------------------------------------------------------
    # Grad (differentiation.py).
    # -------------------------------------------------------------------------
    def grad(self, o):
        #print("\n\nVisiting Grad: " + repr(o))

        # Get expression
        derivative_expr, = o.operands()

        # Get components
        components = self.component()

        en = derivative_expr.rank()
        cn = len(components)
        ffc_assert(o.rank() == cn, "Expecting rank of grad expression to match components length.")

        # Get direction of derivative
        if cn == en+1:
            der = components[en]
            self._components.push(components[:en])
        elif cn == en:
            # This happens in 1D, sligtly messy result of defining grad(f) == f.dx(0)
            der = 0
        else:
            error("Unexpected rank %d and component length %d in grad expression." % (en, cn))

        # Add direction to list of derivatives
        self._derivatives.append(der)

        # Visit children to generate the derivative code.
        code = self.visit(derivative_expr)

        # Remove the direction from list of derivatives
        self._derivatives.pop()
        if cn == en+1:
            self._components.pop()
        return code

    # -------------------------------------------------------------------------
    # Coefficient and Constants (function.py).
    # -------------------------------------------------------------------------
    def coefficient(self, o):
        #print("\nVisiting Coefficient: " + repr(o))

        # Map o to object with proper element and numbering
        o = self._function_replace_map[o]

        # Create aux. info.
        components = self.component()
        derivatives = self.derivatives()

        # Check if function is already in cache
        key = (o, components, derivatives, self.restriction, self.avg)
        function_code = self.function_cache.get(key)

        # FIXME: Why does using a code dict from cache make the expression manipulations blow (MemoryError) up later?
        if function_code is None or self.optimise_parameters["optimisation"]:
            # Get auxiliary variables to generate function
            (component, local_elem, local_comp, local_offset,
             ffc_element, transformation, multiindices) = self._get_auxiliary_variables(o, components, derivatives)

            # Check that we don't take derivatives of QuadratureElements.
            is_quad_element = local_elem.family() == "Quadrature"
            ffc_assert(not (derivatives and is_quad_element), \
                       "Derivatives of Quadrature elements are not supported: " + repr(o))

            # Create code for function and add empty tuple to cache dict.
            function_code = {(): self.create_function(o, derivatives, component,
                                                      local_comp, local_offset, ffc_element, is_quad_element,
                                                      transformation, multiindices, self.tdim, self.gdim, self.avg)}

            self.function_cache[key] = function_code

        return function_code

    def constant(self, o):
        #print("\n\nVisiting Constant: " + repr(o))

        # Map o to object with proper element and numbering
        o = self._function_replace_map[o]

        # Safety checks.
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

    def vector_constant(self, o):
        #print("\n\nVisiting VectorConstant: " + repr(o))

        # Map o to object with proper element and numbering
        o = self._function_replace_map[o]

        # Get the component
        components = self.component()

        # Safety checks.
        ffc_assert(len(components) == 1, "VectorConstant expects 1 component index: " + repr(components))

        # We get one component.
        component = components[0]

        # Handle restriction.
        if self.restriction == "-":
            component += o.shape()[0]

        # Let child class create constant symbol
        coefficient = format["coefficient"](o.count(), component)
        return self._create_symbol(coefficient, CONST)

    def tensor_constant(self, o):
        #print("\n\nVisiting TensorConstant: " + repr(o))

        # Map o to object with proper element and numbering
        o = self._function_replace_map[o]

        # Get the components
        components = self.component()

        # Safety checks.
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
    # SpatialCoordinate (geometry.py).
    # -------------------------------------------------------------------------
    def spatial_coordinate(self, o):
        #print "\n\nVisiting SpatialCoordinate:", repr(o)
        #print "\n\nVisiting SpatialCoordinate:", repr(operands)

        # Get the component.
        components = self.component()
        c, = components

        if self.vertex is not None:
            error("Spatial coordinates (x) not implemented for point measure (dP)") # TODO: Implement this, should be just the point.
        else:
            # Generate the appropriate coordinate and update tables.
            coordinate = format["ip coordinates"](self.points, c)
            self._generate_affine_map()
            return self._create_symbol(coordinate, IP)

    # -------------------------------------------------------------------------
    # Indexed (indexed.py).
    # -------------------------------------------------------------------------
    def indexed(self, o):
        #print("\n\nVisiting Indexed:" + repr(o))

        # Get indexed expression and index, map index to current value
        # and update components
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

    def tan(self, o, *operands):
        #print("\n\nVisiting Tan: " + repr(o) + "with operands: " + "\n".join(map(repr,operands)))
        return self._math_function(operands, format["tan"])

    def cosh(self, o, *operands):
        #print("\n\nVisiting Cosh: " + repr(o) + "with operands: " + "\n".join(map(repr,operands)))
        return self._math_function(operands, format["cosh"])

    def sinh(self, o, *operands):
        #print("\n\nVisiting Sinh: " + repr(o) + "with operands: " + "\n".join(map(repr,operands)))
        return self._math_function(operands, format["sinh"])

    def tanh(self, o, *operands):
        #print("\n\nVisiting Tanh: " + repr(o) + "with operands: " + "\n".join(map(repr,operands)))
        return self._math_function(operands, format["tanh"])

    def acos(self, o, *operands):
        #print("\n\nVisiting Acos: " + repr(o) + "with operands: " + "\n".join(map(repr,operands)))
        return self._math_function(operands, format["acos"])

    def asin(self, o, *operands):
        #print("\n\nVisiting Asin: " + repr(o) + "with operands: " + "\n".join(map(repr,operands)))
        return self._math_function(operands, format["asin"])

    def atan(self, o, *operands):
        #print("\n\nVisiting Atan: " + repr(o) + "with operands: " + "\n".join(map(repr,operands)))
        return self._math_function(operands, format["atan"])

    def atan_2(self, o, *operands):
        #print("\n\nVisiting Atan2: " + repr(o) + "with operands: " + "\n".join(map(repr,operands)))
        self.additional_includes_set.add("#include <cmath>")
        return self._atan_2_function(operands, format["atan_2"])

    def erf(self, o, *operands):
        #print("\n\nVisiting Erf: " + repr(o) + "with operands: " + "\n".join(map(repr,operands)))
        return self._math_function(operands, format["erf"])

    def bessel_i(self, o, *operands):
        #print("\n\nVisiting Bessel_I: " + repr(o) + "with operands: " + "\n".join(map(repr,operands)))
        #self.additional_includes_set.add("#include <tr1/cmath>")
        self.additional_includes_set.add("#include <boost/math/special_functions.hpp>")
        return self._bessel_function(operands, format["bessel_i"])

    def bessel_j(self, o, *operands):
        #print("\n\nVisiting Bessel_J: " + repr(o) + "with operands: " + "\n".join(map(repr,operands)))
        #self.additional_includes_set.add("#include <tr1/cmath>")
        self.additional_includes_set.add("#include <boost/math/special_functions.hpp>")
        return self._bessel_function(operands, format["bessel_j"])

    def bessel_k(self, o, *operands):
        #print("\n\nVisiting Bessel_K: " + repr(o) + "with operands: " + "\n".join(map(repr,operands)))
        #self.additional_includes_set.add("#include <tr1/cmath>")
        self.additional_includes_set.add("#include <boost/math/special_functions.hpp>")
        return self._bessel_function(operands, format["bessel_k"])

    def bessel_y(self, o, *operands):
        #print("\n\nVisiting Bessel_Y: " + repr(o) + "with operands: " + "\n".join(map(repr,operands)))
        #self.additional_includes_set.add("#include <tr1/cmath>")
        self.additional_includes_set.add("#include <boost/math/special_functions.hpp>")
        return self._bessel_function(operands, format["bessel_y"])

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

    def cell_avg(self, o):
        ffc_assert(self.avg is None, "Not expecting nested averages.")

        # Just get the first operand, there should only be one.
        expr, = o.operands()

        # Set average marker, visit operand and reset marker
        self.avg = "cell"
        code = self.visit(expr)
        self.avg = None

        return code

    def facet_avg(self, o):
        ffc_assert(self.avg is None, "Not expecting nested averages.")
        ffc_assert(self.entitytype != "cell", "Cannot take facet_avg in a cell integral.")

        # Just get the first operand, there should only be one.
        expr, = o.operands()

        # Set average marker, visit operand and reset marker
        self.avg = "facet"
        code = self.visit(expr)
        self.avg = None

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
    def generate_terms(self, integrand, domain_type):
        "Generate terms for code generation."
        #print integrand
        #print tree_format(integrand, 0, False)
        # Get terms.
        terms = self.visit(integrand)

        f_nzc = format["nonzero columns"](0).split("0")[0]

        # Loop code and add weight and scale factor to value and sort after
        # loop ranges.
        new_terms = {}
        for key, val in terms.items():
            # If value was zero continue.
            if val is None:
                continue
            # Create data.
            value, ops, sets = self._create_entry_data(val, domain_type)
            # Extract nzc columns if any and add to sets.
            used_nzcs = set([int(k[1].split(f_nzc)[1].split("[")[0]) for k in key if f_nzc in k[1]])
            sets.append(used_nzcs)

            # Create loop information and entry from key info and insert into dict.
            loop, entry = self._create_loop_entry(key, f_nzc)
            if not loop in new_terms:
                sets.append({})
                new_terms[loop] = [sets, [(entry, value, ops)]]
            else:
                for i, s in enumerate(sets):
                    new_terms[loop][0][i].update(s)
                new_terms[loop][1].append((entry, value, ops))

        return new_terms

    def _create_loop_entry(self, key, f_nzc):

        indices = {0: format["first free index"],  1: format["second free index"]}

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
            ffc_assert(key[0] == -2 or key[0] == 0, \
                        "Linear forms must be defined using test functions only: " + repr(key))
            index_j, entry, range_j, space_dim_j = key
            loop = ((indices[index_j], 0, range_j),)
            if range_j == 1 and self.optimise_parameters["ignore ones"] and not (f_nzc in entry):
                loop = ()
        elif len(key) == 2:
            # Extract test and trial loops in correct order and check if for is legal.
            key0, key1 = (0, 0)
            for k in key:
                ffc_assert(k[0] in indices, \
                "Bilinear forms must be defined using test and trial functions (index -2, -1, 0, 1): " + repr(k))
                if k[0] == -2 or k[0] == 0:
                    key0 = k
                else:
                    key1 = k
            index_j, entry_j, range_j, space_dim_j = key0
            index_k, entry_k, range_k, space_dim_k = key1

            loop = []
            if not (range_j == 1 and self.optimise_parameters["ignore ones"]) or f_nzc in entry_j:
                loop.append((indices[index_j], 0, range_j))
            if not (range_k == 1 and self.optimise_parameters["ignore ones"]) or f_nzc in entry_k:
                loop.append((indices[index_k], 0, range_k))
            entry = format["add"]([format["mul"]([entry_j, str(space_dim_k)]), entry_k])
            loop = tuple(loop)
        else:
            error("Only rank 0, 1 and 2 tensors are currently supported: " + repr(key))
        # Generate the code line for the entry.
        # Try to evaluate entry ("3*6 + 2" --> "20").
        try:
            entry = str(eval(entry))
        except:
            pass
        return loop, entry

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

    def _atan_2_function(self, operands, format_function):
        error("This function should be implemented by the child class.")

    def _get_auxiliary_variables(self,
                                 ufl_function,
                                 component,
                                 derivatives):
        "Helper function for both Coefficient and Argument."

        # Get UFL element.
        ufl_element = ufl_function.element()

        # Get subelement and the relative (flattened) component (in case we have mixed elements).
        local_comp, local_elem = ufl_element.extract_component(component)
        ffc_assert(len(local_comp) <= 1, "Assuming there are no tensor-valued basic elements.")
        local_comp = local_comp[0] if local_comp else 0

        # Check that component != not () since the UFL component map will turn
        # it into 0, and () does not mean zeroth component in this context.
        if len(component):
            # Map component using component map from UFL. (TODO: inefficient use of this function)
            comp_map, comp_num = build_component_numbering(ufl_element.value_shape(), ufl_element.symmetry())
            component = comp_map[component]

            # Map physical components into reference components
            component, dummy = transform_component(component, 0, ufl_element)

            # Compute the local offset (needed for non-affine mappings).
            local_offset = component - local_comp
        else:
            # Compute the local offset (needed for non-affine mappings).
            local_offset = 0

        # Create FFC element.
        ffc_element = create_element(ufl_element)

        # Assuming that mappings for all basisfunctions are equal
        # (they should be).
        ffc_sub_element = create_element(local_elem)
        transformation = ffc_sub_element.mapping()[0]

        # Generate FFC multi index for derivatives.
        multiindices = FFCMultiIndex([range(self.tdim)]*len(derivatives)).indices

        #print "in create_auxiliary"
        #print "component = ", component
        return (component, local_elem, local_comp, local_offset, ffc_element, transformation, multiindices)

    def _get_current_entity(self):
        if self.entitytype == "cell":
            # If we add macro cell integration, I guess the 'current cell number' would go here?
            return 0
        elif self.entitytype == "facet":
            # Handle restriction through facet.
            return {"+": self.facet0, "-": self.facet1, None: self.facet0}[self.restriction]
        elif self.entitytype == "vertex":
            return self.vertex
        else:
            error("Unknown entity type %s." % self.entitytype)

    def _create_mapping_basis(self, component, deriv, avg, ufl_argument, ffc_element):
        "Create basis name and mapping from given basis_info."
        ffc_assert(ufl_argument in self._function_replace_values, "Expecting ufl_argument to have been mapped prior to this call.")

        # Get string for integration points.
        f_ip = "0" if (avg or self.points == 1) else format["integration points"]
        generate_psi_name = format["psi name"]

        # Only support test and trial functions.
        indices = { 0: format["first free index"],
                    1: format["second free index"]}

        # Check that we have a basis function.
        ffc_assert(ufl_argument.count() in indices, \
                   "Currently, Argument index must be either 0 or 1: " + repr(ufl_argument))

        # Get element counter and loop index.
        element_counter = self.element_map[1 if avg else self.points][ufl_argument.element()]
        loop_index = indices[ufl_argument.count()]

        # Create basis access, we never need to map the entry in the basis table
        # since we will either loop the entire space dimension or the non-zeros.
        basis_access = format["component"]("", [f_ip, loop_index])

        # Offset element space dimension in case of negative restriction,
        # need to use the complete element for offset in case of mixed element.
        space_dim = ffc_element.space_dimension()
        offset = {"+": "", "-": str(space_dim), None: ""}[self.restriction]

        # If we have a restricted function multiply space_dim by two.
        if self.restriction == "+" or self.restriction == "-":
            space_dim *= 2

        # Get current cell entity, with current restriction considered
        entity = self._get_current_entity()

        name = generate_psi_name(element_counter, self.entitytype, entity, component, deriv, avg)
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

    def _create_function_name(self, component, deriv, avg, is_quad_element, ufl_function, ffc_element):
        ffc_assert(ufl_function in self._function_replace_values, "Expecting ufl_function to have been mapped prior to this call.")

        # Get string for integration points.
        f_ip = "0" if (avg or self.points == 1) else format["integration points"]

        # Get the element counter.
        element_counter = self.element_map[1 if avg else self.points][ufl_function.element()]

        # Get current cell entity, with current restriction considered
        entity = self._get_current_entity()

        # Set to hold used nonzero columns
        used_nzcs = set()

        # Create basis name and map to correct basis and get info.
        generate_psi_name = format["psi name"]
        psi_name = generate_psi_name(element_counter, self.entitytype, entity, component, deriv, avg)
        psi_name, non_zeros, zeros, ones = self.name_map[psi_name]

        # If all basis are zero we just return None.
        if zeros and self.optimise_parameters["ignore zero tables"]:
            return self._format_scalar_value(None)[()]

        # Get the index range of the loop index.
        loop_index_range = shape(self.unique_tables[psi_name])[1]
        if loop_index_range > 1:
            # Pick first free index of secondary type
            # (could use primary indices, but it's better to avoid confusion).
            loop_index = format["free indices"][0]

        # If we have a quadrature element we can use the ip number to look
        # up the value directly. Need to add offset in case of components.
        if is_quad_element:
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
                if non_zeros and f_ip == "0":
                    # If we have non zero column mapping but only one value just pick it.
                    # MSA: This should be an exact refactoring of the previous logic,
                    #      but I'm not sure if these lines were originally intended
                    #      here in the quad_element section, or what this even does:
                    coefficient_access = str(non_zeros[1][0])
                else:
                    coefficient_access = f_ip

        elif non_zeros:
            if loop_index_range == 1:
                # If we have non zero column mapping but only one value just pick it.
                coefficient_access = str(non_zeros[1][0])
            else:
                used_nzcs.add(non_zeros[0])
                coefficient_access = format["component"](format["nonzero columns"](non_zeros[0]), loop_index)

        elif loop_index_range == 1:
            # If the loop index range is one we can look up the first component
            # in the coefficient array.
            coefficient_access = "0"

        else:
            # Or just set default coefficient access.
            coefficient_access = loop_index

        # Offset by element space dimension in case of negative restriction.
        offset = {"+": "", "-": str(ffc_element.space_dimension()), None: ""}[self.restriction]
        if offset:
            coefficient_access = format["add"]([coefficient_access, offset])

        # Try to evaluate coefficient access ("3 + 2" --> "5").
        try:
            coefficient_access = str(eval(coefficient_access))
            C_ACCESS = GEO
        except:
            C_ACCESS = IP
        # Format coefficient access
        coefficient = format["coefficient"](str(ufl_function.count()), coefficient_access)

        # Build and cache some function data only if we need the basis
        # MSA: I don't understand the mix of loop index range check and ones check here, but that's how it was.
        if is_quad_element or (loop_index_range == 1 and ones and self.optimise_parameters["ignore ones"]):
            # If we only have ones or if we have a quadrature element we don't need the basis.
            function_symbol_name = coefficient
            F_ACCESS = C_ACCESS

        else:
            # Add basis name to set of used tables and add matrix access.
            # TODO: We should first add this table if the function is used later
            # in the expressions. If some term is multiplied by zero and it falls
            # away there is no need to compute the function value
            self.used_psi_tables.add(psi_name)

            # Create basis access, we never need to map the entry in the basis
            # table since we will either loop the entire space dimension or the
            # non-zeros.
            basis_index = "0" if loop_index_range == 1 else loop_index
            basis_access = format["component"]("", [f_ip, basis_index])
            basis_name = psi_name + basis_access
            # Try to set access to the outermost possible loop
            if f_ip == "0" and basis_access == "0":
                B_ACCESS = GEO
                F_ACCESS = C_ACCESS
            else:
                B_ACCESS = IP
                F_ACCESS = IP

            # Format expression for function
            function_expr = self._create_product([self._create_symbol(basis_name, B_ACCESS)[()],
                                                  self._create_symbol(coefficient, C_ACCESS)[()]])

            # Check if the expression to compute the function value is already in
            # the dictionary of used function. If not, generate a new name and add.
            data = self.function_data.get(function_expr)
            if data is None:
                function_count = len(self.function_data)
                data = (function_count, loop_index_range,
                        self._count_operations(function_expr),
                        psi_name, used_nzcs, ufl_function.element())
                self.function_data[function_expr] = data
            function_symbol_name = format["function value"](data[0])

        # TODO: This access stuff was changed subtly during my refactoring, the
        # X_ACCESS vars is an attempt at making it right, make sure it is correct now!
        return self._create_symbol(function_symbol_name, F_ACCESS)[()]

    def _generate_affine_map(self):
        """Generate psi table for affine map, used by spatial coordinate to map
        integration point to physical element."""

        # TODO: KBO: Perhaps it is better to create a fiat element and tabulate
        # the values at the integration points?
        f_FEA = format["affine map table"]
        f_ip  = format["integration points"]

        affine_map = {1: lambda x: [1.0 - x[0],               x[0]],
                      2: lambda x: [1.0 - x[0] - x[1],        x[0], x[1]],
                      3: lambda x: [1.0 - x[0] - x[1] - x[2], x[0], x[1], x[2]]}

        num_ip = self.points
        w, points = self.quad_weights[num_ip]

        if self.facet0 is not None:
            points = map_facet_points(points, self.facet0)
            name = f_FEA(num_ip, self.facet0)
        elif self.vertex is not None:
            error("Spatial coordinates (x) not implemented for point measure (dP)") # TODO: Implement this, should be just the point.
            #name = f_FEA(num_ip, self.vertex)
        else:
            name = f_FEA(num_ip, 0)

        if name not in self.unique_tables:
            self.unique_tables[name] = array([affine_map[len(p)](p) for p in points])

        if self.coordinate is None:
            ip = f_ip if num_ip > 1 else 0
            r = None if self.facet1 is None else "+"
            self.coordinate = [name, self.gdim, ip, r]

    # -------------------------------------------------------------------------
    # Helper functions for code_generation()
    # -------------------------------------------------------------------------
    def _count_operations(self, expression):
        error("This function should be implemented by the child class.")

    def _create_entry_data(self, val):
        error("This function should be implemented by the child class.")
