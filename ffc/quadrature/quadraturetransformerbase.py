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
#
# First added:  2009-10-13
# Last changed: 2013-01-08

# Python modules.
from itertools import izip
import time
from numpy import shape, array

# UFL Classes.
from ufl.classes import MultiIndex, FixedIndex, Index
from ufl.common import StackDict, Stack
from ufl.permutation import build_component_numbering

# UFL Algorithms.
from ufl.algorithms import propagate_restrictions, Transformer, tree_format
from ufl.algorithms import strip_variables

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
from symbolics import generate_aux_constants
from symbolics import BASIS, IP, GEO, CONST

class QuadratureTransformerBase(Transformer):
#class QuadratureTransformerBase(ReuseTransformer):
    "Transform UFL representation to quadrature code."

    def __init__(self,
                 psi_tables,
                 quad_weights,
                 geo_dim,
                 top_dim,
                 optimise_parameters):

        Transformer.__init__(self)

        # Save optimise_parameters, weights and fiat_elements_map.
        self.optimise_parameters = optimise_parameters

        # Create containers and variables.
        self.used_psi_tables = set()
        self.psi_tables_map = {}
        self.used_weights = set()
        self.quad_weights = quad_weights
        self.used_nzcs = set()
        self.ip_consts = {}
        self.trans_set = set()
        self.functions = {}
        self.function_count = 0
        self.geo_dim = geo_dim
        self.top_dim = top_dim
        self.points = 0
        self.facet0 = None
        self.facet1 = None
        self.restriction = None
        self.coordinate = None
        self.conditionals = {}
        self.additional_includes_set = set()

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
        self.coordinate = None
        self.conditionals = {}
#        # Reset functions and count everytime we generate a new case of facets.
#        self.functions = {}
#        self.function_count = 0

#        # Reset cache
#        self.argument_cache = {}
#        self.function_cache = {}

    def update_points(self, points):
        self.points = points
        self.coordinate = None
        # Reset functions everytime we move to a new quadrature loop
        self.conditionals = {}
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
        error("All derivatives apart from Grad should have been expanded!!")

    def finite_element_base(self, o, *operands):
        print "\n\nVisiting FiniteElementBase: ", repr(o)
        error("FiniteElements must be member of a Argument or Coefficient!!")

    def form(self, o, *operands):
        print "\n\nVisiting Form: ", repr(o)
        error("The transformer only work on a Form integrand, not the Form itself!!")

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
    def facet_normal(self, o,  *operands):
        print "\n\nVisiting FacetNormal: ", repr(o)
        error("This object should be implemented by the child class.")

    def cell_volume(self, o,  *operands):
        print "\n\nVisiting CellVolume: ", repr(o)
        error("This object should be implemented by the child class.")

    def circumradius(self, o,  *operands):
        print "\n\nVisiting Circumeradius: ", repr(o)
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
        if basis is not None and not self.optimise_parameters["optimisation"]:
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
            ffc_error("Unexpected rank %d and component length %d in grad expression." % (en, cn))

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
        if function_code is not None and not self.optimise_parameters["optimisation"]:
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
    # SpatialCoordinate (geometry.py).
    # -------------------------------------------------------------------------
    def spatial_coordinate(self, o, *operands):
        #print "\n\nVisiting SpatialCoordinate:", repr(o)
        #print "\n\nVisiting SpatialCoordinate:", repr(operands)

        # Get the component.
        components = self.component()

        # Safety checks.
        ffc_assert(not operands, "Didn't expect any operands for spatial_coordinate: " + repr(operands))

        ffc_assert(len(components) == 1,
                   " expects 1 component index: " + repr(components))
        c, = components

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

    def acos(self, o, *operands):
        #print("\n\nVisiting Acos: " + repr(o) + "with operands: " + "\n".join(map(repr,operands)))
        return self._math_function(operands, format["acos"])

    def asin(self, o, *operands):
        #print("\n\nVisiting Asin: " + repr(o) + "with operands: " + "\n".join(map(repr,operands)))
        return self._math_function(operands, format["asin"])

    def atan(self, o, *operands):
        #print("\n\nVisiting Atan: " + repr(o) + "with operands: " + "\n".join(map(repr,operands)))
        return self._math_function(operands, format["atan"])

    def erf(self, o, *operands):
        #print("\n\nVisiting Erf: " + repr(o) + "with operands: " + "\n".join(map(repr,operands)))
        return self._math_function(operands, format["erf"])

    def bessel_i(self, o, *operands):
        #print("\n\nVisiting Bessel_I: " + repr(o) + "with operands: " + "\n".join(map(repr,operands)))
#        self.additional_includes_set.add("#include <tr1/cmath>")
        self.additional_includes_set.add("#include <boost/math/tr1.hpp>")
        return self._bessel_function(operands, format["bessel_i"])

    def bessel_j(self, o, *operands):
        #print("\n\nVisiting Bessel_J: " + repr(o) + "with operands: " + "\n".join(map(repr,operands)))
#        self.additional_includes_set.add("#include <tr1/cmath>")
        self.additional_includes_set.add("#include <boost/math/tr1.hpp>")
        return self._bessel_function(operands, format["bessel_j"])

    def bessel_k(self, o, *operands):
        #print("\n\nVisiting Bessel_K: " + repr(o) + "with operands: " + "\n".join(map(repr,operands)))
#        self.additional_includes_set.add("#include <tr1/cmath>")
        self.additional_includes_set.add("#include <boost/math/tr1.hpp>")
        return self._bessel_function(operands, format["bessel_k"])

    def bessel_y(self, o, *operands):
        #print("\n\nVisiting Bessel_Y: " + repr(o) + "with operands: " + "\n".join(map(repr,operands)))
#        self.additional_includes_set.add("#include <tr1/cmath>")
        self.additional_includes_set.add("#include <boost/math/tr1.hpp>")
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
            value, ops, sets = self._create_entry_data(val)
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

        # TODO: Verify that test and trial functions will ALWAYS be rearranged to 0 and 1.
        indices = {-2: format["first free index"], -1: format["second free index"],
                    0: format["first free index"],  1: format["second free index"]}

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

    def _get_auxiliary_variables(self,
                                 ufl_function,
                                 component,
                                 derivatives):
        "Helper function for both Coefficient and Argument."

        # Get UFL element.
        ufl_element = ufl_function.element()

        # Get local component (in case we have mixed elements).
        local_comp, local_elem = ufl_element.extract_component(component)

        # Check that we don't take derivatives of QuadratureElements.
        quad_element = local_elem.family() == "Quadrature"
        ffc_assert(not (derivatives and quad_element), \
                   "Derivatives of Quadrature elements are not supported: " + repr(ufl_function))

        # Create FFC element.
        ffc_element = create_element(ufl_element)

        # Get relevant sub element and mapping.
        sub_element = create_element(local_elem)

        # Assuming that mappings for all basisfunctions are equal
        # (they should be).
        transformation = sub_element.mapping()[0]

        # Handle tensor elements.
        if len(local_comp) > 1:
            local_comp = local_elem._sub_element_mapping[local_comp]
        elif local_comp:
            local_comp = local_comp[0]
        else:
            local_comp = 0

        # Check that component != not () since the UFL component map will turn
        # it into 0, and () does not mean zeroth component in this context.
        if component != ():
            # Map component using component map from UFL.
            comp_map, comp_num = build_component_numbering(ufl_element.value_shape(), ufl_element.symmetry())
            component = comp_map[component]

        # Map physical components into reference components
        component, dummy = transform_component(component, 0, ufl_element)

        # Compute the local offset (needed for non-affine mappings).
        local_offset = 0
        if component:
            local_offset = component - local_comp

        # Generate FFC multi index for derivatives.
        multiindices = FFCMultiIndex([range(self.top_dim)]*len(derivatives)).indices

        #print "in create_auxiliary"
        #print "component = ", component
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
                data = (self.function_count, loop_index_range, self._count_operations(function_expr),\
                        psi_name, used_nzcs, ufl_function.element())
                self.functions[function_expr] = data
                # Increase count.
                self.function_count += 1
            else:
                data = self.functions[function_expr]
                function_name = self._create_symbol(format["function value"](data[0]), ACCESS)[()]
                # Check just to make sure.
                ffc_assert(data[1] == loop_index_range, "Index ranges does not match." + repr(data[1]) + repr(loop_index_range))
        return function_name

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
        if not self.facet0 is None:
            points = map_facet_points(points, self.facet0)
            name = f_FEA(num_ip, self.facet0)
        else:
            name = f_FEA(num_ip, 0)

        if name not in self.unique_tables:
            vals = []
            for p in points:
                vals.append(affine_map[len(p)](p))
            self.unique_tables[name] = array(vals)
        if self.coordinate is None:
            ip = 0
            r = None
            if num_ip > 1:
                ip = f_ip
            if self.facet1 is not None:
                r = "+"
            self.coordinate = [name, self.geo_dim, ip, r]

    # -------------------------------------------------------------------------
    # Helper functions for code_generation()
    # -------------------------------------------------------------------------
    def _count_operations(self, expression):
        error("This function should be implemented by the child class.")

    def _create_entry_data(self, val):
        error("This function should be implemented by the child class.")
