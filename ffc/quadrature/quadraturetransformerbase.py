# -*- coding: utf-8 -*-
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
# Modified by Martin Sandve Alnæs, 2013
# Modified by Garth N. Wells, 2013
# Modified by Lizao Li, 2015
# Modified by Anders Logg, 2015

# Python modules
import functools
from numpy import shape, array

# UFL Classes
from ufl.classes import FixedIndex, Index
from ufl.utils.stacks import StackDict, Stack
from ufl.permutation import build_component_numbering
from ufl import custom_integral_types

# UFL Algorithms
from ufl.algorithms import Transformer

# FFC modules.
from ffc.log import ffc_assert, error
from ffc.fiatinterface import MixedElement, create_element, map_facet_points
from ffc.quadrature.cpp import format

# FFC utils
from ffc.utils import listcopy
from ffc.representationutils import transform_component

# Utility and optimisation functions for quadraturegenerator.
from ffc.quadrature.quadratureutils import create_psi_tables
from ffc.quadrature.symbolics import BASIS, IP, GEO


class FFCMultiIndex:
    """A MultiIndex represents a list of indices and holds the following
    data:

        rank    - rank of multiindex
        dims    - a list of dimensions
        indices - a list of all possible multiindex values

    """

    def __init__(self, dims):
        "Create multiindex from given list of ranges"

        def outer_join(a, b):
            """Let a be a list of lists and b a list. We append each element of b
            to each list in a and return the resulting list of lists.

            """
            outer = []
            for i in range(len(a)):
                for j in range(len(b)):
                    outer += [a[i] + [b[j]]]
            return outer

        def build_indices(dims):
            "Create a list of all index combinations."
            if not dims:
                return [[]]
            ranges = listcopy(dims)
            return functools.reduce(outer_join, ranges, [[]])

        self.rank = len(dims)
        self.dims = [len(dim) for dim in dims]
        self.indices = build_indices(dims)
        return

    def __str__(self):
        "Return informal string representation (pretty-print)."
        return "rank = %d dims = %s indices = %s" % (self.rank, str(self.dims), str(self.indices))


class QuadratureTransformerBase(Transformer):

    "Transform UFL representation to quadrature code."

    def __init__(self,
                 psi_tables,
                 quad_weights,
                 gdim,
                 tdim,
                 entity_type,
                 function_replace_map,
                 optimise_parameters):

        Transformer.__init__(self)

        # Save optimise_parameters, weights and fiat_elements_map.
        self.optimise_parameters = optimise_parameters

        # Map from original functions with possibly incomplete
        # elements to functions with properly completed elements
        self._function_replace_map = function_replace_map
        self._function_replace_values = set(function_replace_map.values())  # For assertions

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
        self.entity_type = entity_type
        self.points = 0
        self.facet0 = None
        self.facet1 = None
        self.vertex = None
        self.restriction = None
        self.avg = None
        self.coordinate = None
        self.conditionals = {}
        self.additional_includes_set = set()

        # Stacks.
        self._derivatives = []
        self._index2value = StackDict()
        self._components = Stack()

        self.element_map, self.name_map, self.unique_tables =\
            create_psi_tables(psi_tables,
                              self.optimise_parameters["eliminate zeros"],
                              self.entity_type)

        # Cache.
        self.argument_cache = {}
        self.function_cache = {}

    def update_cell(self):
        ffc_assert(self.entity_type == "cell",
                   "Not expecting update_cell on a %s." % self.entity_type)
        self.facet0 = None
        self.facet1 = None
        self.vertex = None
        self.coordinate = None
        self.conditionals = {}

    def update_facets(self, facet0, facet1):
        ffc_assert(self.entity_type == "facet",
                   "Not expecting update_facet on a %s." % self.entity_type)
        self.facet0 = facet0
        self.facet1 = facet1
        self.vertex = None
        self.coordinate = None
        self.conditionals = {}

    def update_vertex(self, vertex):
        ffc_assert(self.entity_type == "vertex",
                   "Not expecting update_vertex on a %s." % self.entity_type)
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
        print("\n\n **** Displaying QuadratureTransformer ****")
        print("\nQuadratureTransformer, element_map:\n", self.element_map)
        print("\nQuadratureTransformer, name_map:\n", self.name_map)
        print("\nQuadratureTransformer, unique_tables:\n", self.unique_tables)
        print("\nQuadratureTransformer, used_psi_tables:\n",
              self.used_psi_tables)
        print("\nQuadratureTransformer, psi_tables_map:\n", self.psi_tables_map)
        print("\nQuadratureTransformer, used_weights:\n", self.used_weights)

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
        print("\n\nVisiting basic Expr:", repr(o), "with operands:")
        error("This expression is not handled: " + repr(o))

    # Nothing in terminal.py is handled. Can only handle children of
    # these clases.
    def terminal(self, o):
        print("\n\nVisiting basic Terminal:", repr(o), "with operands:")
        error("This terminal is not handled: " + repr(o))

    # -------------------------------------------------------------------------
    # Things which should not be here (after expansion etc.) from:
    # algebra.py, differentiation.py, finiteelement.py, form.py,
    # geometry.py, indexing.py, integral.py, tensoralgebra.py,
    # variable.py.
    # -------------------------------------------------------------------------
    def derivative(self, o, *operands):
        print("\n\nVisiting Derivative: ", repr(o))
        error("All derivatives apart from Grad should have been expanded!!")

    def compound_tensor_operator(self, o):
        print("\n\nVisiting CompoundTensorOperator: ", repr(o))
        error("CompoundTensorOperator should have been expanded.")

    def label(self, o):
        print("\n\nVisiting Label: ", repr(o))
        error("What is a Lable doing in the integrand?")

    # -------------------------------------------------------------------------
    # Things which are not supported yet, from:
    # condition.py, constantvalue.py, function.py, geometry.py, lifting.py,
    # mathfunctions.py, restriction.py
    # -------------------------------------------------------------------------
    def condition(self, o):
        print("\n\nVisiting Condition:", repr(o))
        error("This type of Condition is not supported (yet).")

    def constant_value(self, o):
        print("\n\nVisiting ConstantValue:", repr(o))
        error("This type of ConstantValue is not supported (yet).")

    def geometric_quantity(self, o):
        print("\n\nVisiting GeometricQuantity:", repr(o))
        error("This type of GeometricQuantity is not supported (yet).")

    def math_function(self, o):
        print("\n\nVisiting MathFunction:", repr(o))
        error("This MathFunction is not supported (yet).")

    def atan_2_function(self, o):
        print("\n\nVisiting Atan2Function:", repr(o))
        error("Atan2Function is not implemented (yet).")

    def bessel_function(self, o):
        print("\n\nVisiting BesselFunction:", repr(o))
        error("BesselFunction is not implemented (yet).")

    def restricted(self, o):
        print("\n\nVisiting Restricted:", repr(o))
        error("This type of Restricted is not supported (only positive and negative are currently supported).")

    # -------------------------------------------------------------------------
    # Handlers that should be implemented by child classes.
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    # AlgebraOperators (algebra.py).
    # -------------------------------------------------------------------------
    def sum(self, o, *operands):
        print("\n\nVisiting Sum: ", repr(o))
        error("This object should be implemented by the child class.")

    def product(self, o, *operands):
        print("\n\nVisiting Product: ", repr(o))
        error("This object should be implemented by the child class.")

    def division(self, o, *operands):
        print("\n\nVisiting Division: ", repr(o))
        error("This object should be implemented by the child class.")

    def power(self, o):
        print("\n\nVisiting Power: ", repr(o))
        error("This object should be implemented by the child class.")

    def abs(self, o, *operands):
        print("\n\nVisiting Abs: ", repr(o))
        error("This object should be implemented by the child class.")

    # -------------------------------------------------------------------------
    # FacetNormal, CellVolume, Circumradius (geometry.py).
    # -------------------------------------------------------------------------
    def cell_coordinate(self, o):
        error("This object should be implemented by the child class.")

    def facet_coordinate(self, o):
        error("This object should be implemented by the child class.")

    def cell_origin(self, o):
        error("This object should be implemented by the child class.")

    def facet_origin(self, o):
        error("This object should be implemented by the child class.")

    def cell_facet_origin(self, o):
        error("This object should be implemented by the child class.")

    def jacobian(self, o):
        error("This object should be implemented by the child class.")

    def jacobian_determinant(self, o):
        error("This object should be implemented by the child class.")

    def jacobian_inverse(self, o):
        error("This object should be implemented by the child class.")

    def facet_jacobian(self, o):
        error("This object should be implemented by the child class.")

    def facet_jacobian_determinant(self, o):
        error("This object should be implemented by the child class.")

    def facet_jacobian_inverse(self, o):
        error("This object should be implemented by the child class.")

    def cell_facet_jacobian(self, o):
        error("This object should be implemented by the child class.")

    def cell_facet_jacobian_determinant(self, o):
        error("This object should be implemented by the child class.")

    def cell_facet_jacobian_inverse(self, o):
        error("This object should be implemented by the child class.")

    def facet_normal(self, o):
        error("This object should be implemented by the child class.")

    def cell_normal(self, o):
        error("This object should be implemented by the child class.")

    def cell_volume(self, o):
        error("This object should be implemented by the child class.")

    def circumradius(self, o):
        error("This object should be implemented by the child class.")

    def facet_area(self, o):
        error("This object should be implemented by the child class.")

    def min_facet_edge_length(self, o):
        error("This object should be implemented by the child class.")

    def max_facet_edge_length(self, o):
        error("This object should be implemented by the child class.")

    def cell_orientation(self, o):
        error("This object should be implemented by the child class.")

    def quadrature_weight(self, o):
        error("This object should be implemented by the child class.")

    # -------------------------------------------------------------------------
    # Things that can be handled by the base class.
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    # Argument (basisfunction.py).
    # -------------------------------------------------------------------------
    def argument(self, o):
        # print("\nVisiting Argument:" + repr(o))

        # Create aux. info.
        components = self.component()
        derivatives = self.derivatives()

        # Check if basis is already in cache
        key = (o, components, derivatives, self.restriction, self.avg)
        basis = self.argument_cache.get(key, None)

        tdim = self.tdim

        # FIXME: Why does using a code dict from cache make the
        # expression manipulations blow (MemoryError) up later?
        if basis is None or self.optimise_parameters["optimisation"]:
            # Get auxiliary variables to generate basis
            (component, local_elem, local_comp, local_offset, ffc_element,
             transformation,
             multiindices) = self._get_auxiliary_variables(o, components,
                                                           derivatives)

            # Create mapping and code for basis function and add to
            # dict.
            basis = self.create_argument(o, derivatives, component, local_comp,
                                         local_offset, ffc_element,
                                         transformation, multiindices,
                                         tdim, self.gdim, self.avg)
            self.argument_cache[key] = basis

        return basis

    # -------------------------------------------------------------------------
    # Constant values (constantvalue.py).
    # -------------------------------------------------------------------------
    def identity(self, o):

        # Get components
        i, j = self.component()

        # Only return a value if i==j
        if i == j:
            return self._format_scalar_value(1.0)
        else:
            return self._format_scalar_value(None)

    def scalar_value(self, o):
        "ScalarValue covers IntValue and FloatValue"
        return self._format_scalar_value(o.value())

    def zero(self, o):
        return self._format_scalar_value(None)

    # -------------------------------------------------------------------------
    # Grad (differentiation.py).
    # -------------------------------------------------------------------------
    def grad(self, o):

        # Get expression
        derivative_expr, = o.ufl_operands

        # Get components
        components = self.component()

        en = len(derivative_expr.ufl_shape)
        cn = len(components)
        ffc_assert(len(o.ufl_shape) == cn,
                   "Expecting rank of grad expression to match components length.")

        # Get direction of derivative
        if cn == en + 1:
            der = components[en]
            self._components.push(components[:en])
        elif cn == en:
            # This happens in 1D, sligtly messy result of defining
            # grad(f) == f.dx(0)
            der = 0
        else:
            error("Unexpected rank %d and component length %d in grad expression." % (en, cn))

        # Add direction to list of derivatives
        self._derivatives.append(der)

        # Visit children to generate the derivative code.
        code = self.visit(derivative_expr)

        # Remove the direction from list of derivatives
        self._derivatives.pop()
        if cn == en + 1:
            self._components.pop()
        return code

    # -------------------------------------------------------------------------
    # Coefficient and Constants (function.py).
    # -------------------------------------------------------------------------
    def coefficient(self, o):
        # print("\nVisiting Coefficient: " + repr(o))

        # Map o to object with proper element and count
        o = self._function_replace_map[o]

        # Create aux. info.
        components = self.component()
        derivatives = self.derivatives()

        # Check if function is already in cache
        key = (o, components, derivatives, self.restriction, self.avg)
        function_code = self.function_cache.get(key)

        # FIXME: Why does using a code dict from cache make the
        # expression manipulations blow (MemoryError) up later?
        if function_code is None or self.optimise_parameters["optimisation"]:
            # Get auxiliary variables to generate function
            (component, local_elem, local_comp, local_offset,
             ffc_element, transformation,
             multiindices) = self._get_auxiliary_variables(o, components,
                                                           derivatives)

            # Check that we don't take derivatives of
            # QuadratureElements.
            is_quad_element = local_elem.family() == "Quadrature"
            ffc_assert(not (derivatives and is_quad_element),
                       "Derivatives of Quadrature elements are not supported: " + repr(o))

            tdim = self.tdim

            # Create code for function and add empty tuple to cache
            # dict.
            function_code = {(): self.create_function(o, derivatives, component,
                                                      local_comp, local_offset,
                                                      ffc_element,
                                                      is_quad_element,
                                                      transformation,
                                                      multiindices, tdim,
                                                      self.gdim, self.avg)}

            self.function_cache[key] = function_code

        return function_code

    # -------------------------------------------------------------------------
    # SpatialCoordinate (geometry.py).
    # -------------------------------------------------------------------------
    def spatial_coordinate(self, o):

        # Get the component.
        components = self.component()
        c, = components

        if self.vertex is not None:
            error("Spatial coordinates (x) not implemented for point measure (dP)")  # TODO: Implement this, should be just the point.
        elif self.points is None:
            gdim, = o.ufl_shape
            coordinate = "quadrature_points[ip*%d + %d]" % (gdim, c)
            return self._create_symbol(coordinate, IP)
        else:
            # Generate the appropriate coordinate and update tables.
            coordinate = format["ip coordinates"](self.points, c)
            self._generate_affine_map()
            return self._create_symbol(coordinate, IP)

    # -------------------------------------------------------------------------
    # Indexed (indexed.py).
    # -------------------------------------------------------------------------
    def indexed(self, o):

        # Get indexed expression and index, map index to current value
        # and update components
        indexed_expr, index = o.ufl_operands
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

        # Get expression and index that we're summing over
        summand, multiindex = o.ufl_operands
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
        self.additional_includes_set.add("#include <cmath>")
        return self._math_function(operands, format["sqrt"])

    def exp(self, o, *operands):
        self.additional_includes_set.add("#include <cmath>")
        return self._math_function(operands, format["exp"])

    def ln(self, o, *operands):
        self.additional_includes_set.add("#include <cmath>")
        return self._math_function(operands, format["ln"])

    def cos(self, o, *operands):
        self.additional_includes_set.add("#include <cmath>")
        return self._math_function(operands, format["cos"])

    def sin(self, o, *operands):
        self.additional_includes_set.add("#include <cmath>")
        return self._math_function(operands, format["sin"])

    def tan(self, o, *operands):
        self.additional_includes_set.add("#include <cmath>")
        return self._math_function(operands, format["tan"])

    def cosh(self, o, *operands):
        self.additional_includes_set.add("#include <cmath>")
        return self._math_function(operands, format["cosh"])

    def sinh(self, o, *operands):
        self.additional_includes_set.add("#include <cmath>")
        return self._math_function(operands, format["sinh"])

    def tanh(self, o, *operands):
        self.additional_includes_set.add("#include <cmath>")
        return self._math_function(operands, format["tanh"])

    def acos(self, o, *operands):
        self.additional_includes_set.add("#include <cmath>")
        return self._math_function(operands, format["acos"])

    def asin(self, o, *operands):
        self.additional_includes_set.add("#include <cmath>")
        return self._math_function(operands, format["asin"])

    def atan(self, o, *operands):
        self.additional_includes_set.add("#include <cmath>")
        return self._math_function(operands, format["atan"])

    def atan_2(self, o, *operands):
        self.additional_includes_set.add("#include <cmath>")
        return self._atan_2_function(operands, format["atan_2"])

    def erf(self, o, *operands):
        self.additional_includes_set.add("#include <cmath>")
        return self._math_function(operands, format["erf"])

    def bessel_i(self, o, *operands):
        self.additional_includes_set.add("#include <boost/math/special_functions.hpp>")
        return self._bessel_function(operands, format["bessel_i"])

    def bessel_j(self, o, *operands):
        self.additional_includes_set.add("#include <boost/math/special_functions.hpp>")
        return self._bessel_function(operands, format["bessel_j"])

    def bessel_k(self, o, *operands):
        self.additional_includes_set.add("#include <boost/math/special_functions.hpp>")
        return self._bessel_function(operands, format["bessel_k"])

    def bessel_y(self, o, *operands):
        self.additional_includes_set.add("#include <boost/math/special_functions.hpp>")
        return self._bessel_function(operands, format["bessel_y"])

    # -------------------------------------------------------------------------
    # PositiveRestricted and NegativeRestricted (restriction.py).
    # -------------------------------------------------------------------------
    def positive_restricted(self, o):

        # Just get the first operand, there should only be one.
        restricted_expr = o.ufl_operands
        ffc_assert(len(restricted_expr) == 1,
                   "Only expected one operand for restriction: " + repr(restricted_expr))
        ffc_assert(self.restriction is None,
                   "Expression is restricted twice: " + repr(restricted_expr))

        # Set restriction, visit operand and reset restriction
        self.restriction = "+"
        code = self.visit(restricted_expr[0])
        self.restriction = None

        return code

    def negative_restricted(self, o):

        # Just get the first operand, there should only be one.
        restricted_expr = o.ufl_operands
        ffc_assert(len(restricted_expr) == 1,
                   "Only expected one operand for restriction: " + repr(restricted_expr))
        ffc_assert(self.restriction is None,
                   "Expression is restricted twice: " + repr(restricted_expr))

        # Set restriction, visit operand and reset restriction
        self.restriction = "-"
        code = self.visit(restricted_expr[0])
        self.restriction = None

        return code

    def cell_avg(self, o):
        ffc_assert(self.avg is None, "Not expecting nested averages.")

        # Just get the first operand, there should only be one.
        expr, = o.ufl_operands

        # Set average marker, visit operand and reset marker
        self.avg = "cell"
        code = self.visit(expr)
        self.avg = None

        return code

    def facet_avg(self, o):
        ffc_assert(self.avg is None, "Not expecting nested averages.")
        ffc_assert(self.entity_type != "cell",
                   "Cannot take facet_avg in a cell integral.")

        # Just get the first operand, there should only be one.
        expr, = o.ufl_operands

        # Set average marker, visit operand and reset marker
        self.avg = "facet"
        code = self.visit(expr)
        self.avg = None

        return code

    # -------------------------------------------------------------------------
    # ComponentTensor (tensors.py).
    # -------------------------------------------------------------------------
    def component_tensor(self, o):

        # Get expression and indices
        component_expr, indices = o.ufl_operands

        # Get current component(s)
        components = self.component()

        ffc_assert(len(components) == len(indices),
                   "The number of known components must be equal to the number of components of the ComponentTensor for this to work.")

        # Update the index dict (map index values of current known
        # indices to those of the component tensor)
        for i, v in zip(indices._indices, components):
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

        # Get the component
        component = self.component()

        # Extract first and the rest of the components
        c0, c1 = component[0], component[1:]

        # Get first operand
        op = o.ufl_operands[c0]

        # Evaluate subtensor with this subcomponent
        self._components.push(c1)
        code = self.visit(op)
        self._components.pop()

        return code

    # -------------------------------------------------------------------------
    # Variable (variable.py).
    # -------------------------------------------------------------------------
    def variable(self, o):
        return self.visit(o.expression())

    # -------------------------------------------------------------------------
    # Generate terms for representation.
    # -------------------------------------------------------------------------
    def generate_terms(self, integrand, integral_type):
        "Generate terms for code generation."

        # Set domain type
        self.integral_type = integral_type

        # Get terms
        terms = self.visit(integrand)

        # Get formatting
        f_nzc = format["nonzero columns"](0).split("0")[0]

        # Loop code and add weight and scale factor to value and sort
        # after loop ranges.
        new_terms = {}
        for key, val in sorted(terms.items()):
            # If value was zero continue.
            if val is None:
                continue
            # Create data.
            value, ops, sets = self._create_entry_data(val, integral_type)
            # Extract nzc columns if any and add to sets.
            used_nzcs = set([int(k[1].split(f_nzc)[1].split("[")[0]) for k in key if f_nzc in k[1]])
            sets.append(used_nzcs)

            # Create loop information and entry from key info and
            # insert into dict.
            loop, entry = self._create_loop_entry(key, f_nzc)
            if loop not in new_terms:
                sets.append({})
                new_terms[loop] = [sets, [(entry, value, ops)]]
            else:
                for i, s in enumerate(sets):
                    new_terms[loop][0][i].update(s)
                new_terms[loop][1].append((entry, value, ops))

        return new_terms

    def _create_loop_entry(self, key, f_nzc):

        indices = {0: format["first free index"],
                   1: format["second free index"]}

        # Create appropriate entries.
        # FIXME: We only support rank 0, 1 and 2.
        entry = ""
        loop = ()
        if len(key) == 0:
            entry = "0"
        elif len(key) == 1:
            key = key[0]
            # Checking if the basis was a test function.
            # TODO: Make sure test function indices are always
            # rearranged to 0.
            ffc_assert(key[0] == -2 or key[0] == 0,
                       "Linear forms must be defined using test functions only: " + repr(key))
            index_j, entry, range_j, space_dim_j = key
            loop = ((indices[index_j], 0, range_j),)
            if range_j == 1 and self.optimise_parameters["ignore ones"] and not (f_nzc in entry):
                loop = ()
        elif len(key) == 2:
            # Extract test and trial loops in correct order and check
            # if for is legal.
            key0, key1 = (0, 0)
            for k in key:
                ffc_assert(k[0] in indices,
                           "Bilinear forms must be defined using test and trial functions (index -2, -1, 0, 1): " + repr(k))
                if k[0] == -2 or k[0] == 0:
                    key0 = k
                else:
                    key1 = k
            index_j, entry_j, range_j, space_dim_j = key0
            index_k, entry_k, range_k, space_dim_k = key1

            loop = []
            if not (range_j == 1 and
                    self.optimise_parameters["ignore ones"]) or f_nzc in entry_j:
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
        except Exception:
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

    def _get_auxiliary_variables(self, ufl_function, component, derivatives):
        "Helper function for both Coefficient and Argument."

        # Get UFL element.
        ufl_element = ufl_function.ufl_element()

        # Get subelement and the relative (flattened) component (in
        # case we have mixed elements).
        local_comp, local_elem = ufl_element.extract_component(component)

        # For basic tensor elements, local_comp should be flattened
        if len(local_comp) and len(local_elem.value_shape()) > 0:
            # Map component using component map from UFL. (TODO:
            # inefficient use of this function)
            comp_map, _ = build_component_numbering(local_elem.value_shape(),
                                                    local_elem.symmetry())
            local_comp = comp_map[local_comp]

        # Set local_comp to 0 if it is ()
        if not local_comp:
            local_comp = 0

        # Check that component != not () since the UFL component map
        # will turn it into 0, and () does not mean zeroth component
        # in this context.
        if len(component):
            # Map component using component map from UFL. (TODO:
            # inefficient use of this function)
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
        ffc_sub_element = create_element(local_elem)
        transformation = ffc_sub_element.mapping()[0]
        ffc_assert(all(transformation == mapping for mapping in ffc_sub_element.mapping()),
                   "Assuming subelement mappings are equal but they differ.")

        # Generate FFC multi index for derivatives.
        tdim = self.tdim
        multiindices = FFCMultiIndex([list(range(tdim))] * len(derivatives)).indices

        return (component, local_elem, local_comp, local_offset, ffc_element,
                transformation, multiindices)

    def _get_current_entity(self):
        if self.entity_type == "cell":
            # If we add macro cell integration, I guess the 'current
            # cell number' would go here?
            return 0
        elif self.entity_type == "facet":
            # Handle restriction through facet.
            return {"+": self.facet0, "-": self.facet1,
                    None: self.facet0}[self.restriction]
        elif self.entity_type == "vertex":
            return self.vertex
        else:
            error("Unknown entity type %s." % self.entity_type)

    def _create_mapping_basis(self, component, deriv, avg, ufl_argument,
                              ffc_element):
        "Create basis name and mapping from given basis_info."

        # Get string for integration points.
        f_ip = "0" if (avg or self.points == 1) else format["integration points"]
        generate_psi_name = format["psi name"]

        # Only support test and trial functions.
        indices = {0: format["first free index"],
                   1: format["second free index"]}

        # Check that we have a basis function.
        ffc_assert(ufl_argument.number() in indices,
                   "Currently, Argument number must be either 0 or 1: " + repr(ufl_argument))
        ffc_assert(ufl_argument.part() is None,
                   "Currently, Argument part is not supporte: " + repr(ufl_argument))

        # Get element counter and loop index.
        element_counter = self.element_map[1 if avg else self.points][ufl_argument.ufl_element()]
        loop_index = indices[ufl_argument.number()]

        # Offset element space dimension in case of negative
        # restriction, need to use the complete element for offset in
        # case of mixed element.
        space_dim = ffc_element.space_dimension()
        offset = {"+": "", "-": str(space_dim), None: ""}[self.restriction]

        # If we have a restricted function multiply space_dim by two.
        if self.restriction in ("+", "-"):
            space_dim *= 2

        # Create basis access, we never need to map the entry in the
        # basis table since we will either loop the entire space
        # dimension or the non-zeros.
        if self.restriction in ("+", "-") and self.integral_type in custom_integral_types and offset != "":
            # Special case access for custom integrals (all basis
            # functions stored in flattened array)
            basis_access = format["component"]("", [f_ip, format["add"]([loop_index, offset])])
        else:
            # Normal basis function access
            basis_access = format["component"]("", [f_ip, loop_index])

        # Get current cell entity, with current restriction considered
        entity = self._get_current_entity()
        name = generate_psi_name(element_counter, self.entity_type, entity,
                                 component, deriv, avg)
        name, non_zeros, zeros, ones = self.name_map[name]
        loop_index_range = shape(self.unique_tables[name])[1]

        # If domain type is custom, then special-case set loop index
        # range since table is empty
        if self.integral_type in custom_integral_types:
            loop_index_range = ffc_element.space_dimension()  # different from `space_dimension`...

        basis = ""
        # Ignore zeros if applicable
        if zeros and (self.optimise_parameters["ignore zero tables"] or self.optimise_parameters["remove zero terms"]):
            basis = self._format_scalar_value(None)[()]
        # If the loop index range is one we can look up the first
        # component in the psi array. If we only have ones we don't
        # need the basis.
        elif self.optimise_parameters["ignore ones"] and loop_index_range == 1 and ones:
            loop_index = "0"
            basis = self._format_scalar_value(1.0)[()]
        else:
            # Add basis name to the psi tables map for later use.
            basis = self._create_symbol(name + basis_access, BASIS)[()]
            self.psi_tables_map[basis] = name

        # Create the correct mapping of the basis function into the
        # local element tensor.
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
        except Exception:
            pass

        # Create mapping (index, map, loop_range, space_dim).
        # Example dx and ds: (0, j, 3, 3)
        # Example dS: (0, (j + 3), 3, 6), 6=2*space_dim
        # Example dS optimised: (0, (nz2[j] + 3), 2, 6), 6=2*space_dim
        mapping = ((ufl_argument.number(), basis_map, loop_index_range,
                    space_dim),)

        return (mapping, basis)

    def _create_function_name(self, component, deriv, avg, is_quad_element,
                              ufl_function, ffc_element):
        ffc_assert(ufl_function in self._function_replace_values,
                   "Expecting ufl_function to have been mapped prior to this call.")

        # Get string for integration points.
        f_ip = "0" if (avg or self.points == 1) else format["integration points"]

        # Get the element counter.
        element_counter = self.element_map[1 if avg else self.points][ufl_function.ufl_element()]

        # Get current cell entity, with current restriction considered
        entity = self._get_current_entity()

        # Set to hold used nonzero columns
        used_nzcs = set()

        # Create basis name and map to correct basis and get info.
        generate_psi_name = format["psi name"]
        psi_name = generate_psi_name(element_counter, self.entity_type, entity,
                                     component, deriv, avg)
        psi_name, non_zeros, zeros, ones = self.name_map[psi_name]

        # If all basis are zero we just return None.
        if zeros and self.optimise_parameters["ignore zero tables"]:
            return self._format_scalar_value(None)[()]

        # Get the index range of the loop index.
        loop_index_range = shape(self.unique_tables[psi_name])[1]

        # If domain type is custom, then special-case set loop index
        # range since table is empty
        if self.integral_type in custom_integral_types:
            loop_index_range = ffc_element.space_dimension()

        # Create loop index
        if loop_index_range > 1:
            # Pick first free index of secondary type (could use
            # primary indices, but it's better to avoid confusion).
            loop_index = format["free indices"][0]

        # If we have a quadrature element we can use the ip number to look
        # up the value directly. Need to add offset in case of components.
        if is_quad_element:
            quad_offset = 0
            if component:
                # FIXME: Should we add a member function elements() to
                # FiniteElement?
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
                    # If we have non zero column mapping but only one
                    # value just pick it.
                    # MSA: This should be an exact refactoring of the
                    #      previous logic, but I'm not sure if these
                    #      lines were originally intended here in the
                    #      quad_element section, or what this even
                    #      does:
                    coefficient_access = str(non_zeros[1][0])
                else:
                    coefficient_access = f_ip

        elif non_zeros:
            if loop_index_range == 1:
                # If we have non zero column mapping but only one
                # value just pick it.
                coefficient_access = str(non_zeros[1][0])
            else:
                used_nzcs.add(non_zeros[0])
                coefficient_access = format["component"](format["nonzero columns"](non_zeros[0]), loop_index)

        elif loop_index_range == 1:
            # If the loop index range is one we can look up the first
            # component in the coefficient array.
            coefficient_access = "0"

        else:
            # Or just set default coefficient access.
            coefficient_access = loop_index

        # Offset by element space dimension in case of negative
        # restriction.
        offset = {"+": "", "-": str(ffc_element.space_dimension()), None: ""}[self.restriction]
        if offset:
            coefficient_access = format["add"]([coefficient_access, offset])

        # Try to evaluate coefficient access ("3 + 2" --> "5").
        try:
            coefficient_access = str(eval(coefficient_access))
            C_ACCESS = GEO
        except Exception:
            C_ACCESS = IP
        # Format coefficient access
        coefficient = format["coefficient"](str(ufl_function.count()),
                                            coefficient_access)

        # Build and cache some function data only if we need the basis
        # MSA: I don't understand the mix of loop index range check
        # and ones check here, but that's how it was.
        if is_quad_element or (loop_index_range == 1 and ones and self.optimise_parameters["ignore ones"]):
            # If we only have ones or if we have a quadrature element
            # we don't need the basis.
            function_symbol_name = coefficient
            F_ACCESS = C_ACCESS

        else:
            # Add basis name to set of used tables and add matrix
            # access.
            # TODO: We should first add this table if the function is
            # used later in the expressions. If some term is
            # multiplied by zero and it falls away there is no need to
            # compute the function value
            self.used_psi_tables.add(psi_name)

            # Create basis access, we never need to map the entry in
            # the basis table since we will either loop the entire
            # space dimension or the non-zeros.
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

            # Check if the expression to compute the function value is
            # already in the dictionary of used function. If not,
            # generate a new name and add.
            data = self.function_data.get(function_expr)
            if data is None:
                function_count = len(self.function_data)
                data = (function_count, loop_index_range,
                        self._count_operations(function_expr),
                        psi_name, used_nzcs, ufl_function.ufl_element())
                self.function_data[function_expr] = data
            function_symbol_name = format["function value"](data[0])

        # TODO: This access stuff was changed subtly during my
        # refactoring, the
        # X_ACCESS vars is an attempt at making it right, make sure it
        # is correct now!
        return self._create_symbol(function_symbol_name, F_ACCESS)[()]

    def _generate_affine_map(self):
        """Generate psi table for affine map, used by spatial coordinate to
        map integration point to physical element.

        """

        # TODO: KBO: Perhaps it is better to create a fiat element and
        # tabulate the values at the integration points?
        f_FEA = format["affine map table"]
        f_ip = format["integration points"]

        affine_map = {1: lambda x: [1.0 - x[0], x[0]],
                      2: lambda x: [1.0 - x[0] - x[1], x[0], x[1]],
                      3: lambda x: [1.0 - x[0] - x[1] - x[2], x[0], x[1], x[2]]}

        num_ip = self.points
        w, points = self.quad_weights[num_ip]

        if self.facet0 is not None:
            # Extract the geometric dimension of the points we want to map
            dim = len(points[0]) + 1

            dim2cellname = ["interval", "triangle", "tetrahedron"]
            cellname = dim2cellname[dim-1]

            points = map_facet_points(points, self.facet0, cellname)
            name = f_FEA(num_ip, self.facet0)
        elif self.vertex is not None:
            error("Spatial coordinates (x) not implemented for point measure (dP)")  # TODO: Implement this, should be just the point.
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
