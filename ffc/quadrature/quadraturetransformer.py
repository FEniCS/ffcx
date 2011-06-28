"QuadratureTransformer for quadrature code generation to translate UFL expressions."

# Copyright (C) 2009-2011 Kristian B. Oelgaard
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
# Modified by Peter Brune, 2009
# Modified by Anders Logg, 2009
#
# First added:  2009-02-09
# Last changed: 2011-06-28

# Python modules.
from numpy import shape

# UFL common.
from ufl.common import product, StackDict, Stack

# UFL Classes.
from ufl.classes import FixedIndex
from ufl.classes import IntValue
from ufl.classes import FloatValue
from ufl.classes import Coefficient
from ufl.expr import Operator

# UFL Algorithms.
from ufl.algorithms.printing import tree_format

# FFC modules.
from ffc.log import info, debug, error, ffc_assert
from ffc.cpp import format

# Utility and optimisation functions for quadraturegenerator.
from quadraturetransformerbase import QuadratureTransformerBase
from quadratureutils import create_permutations
from reduce_operations import operation_count
from symbolics import IP

class QuadratureTransformer(QuadratureTransformerBase):
    "Transform UFL representation to quadrature code."

    def __init__(self, *args):

        # Initialise base class.
        QuadratureTransformerBase.__init__(self, *args)

    # -------------------------------------------------------------------------
    # Start handling UFL classes.
    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------
    # AlgebraOperators (algebra.py).
    # -------------------------------------------------------------------------
    def sum(self, o, *operands):
        #print("Visiting Sum: " + "\noperands: \n" + "\n".join(map(repr, operands)))

        # Prefetch formats to speed up code generation.
        f_group  = format["grouping"]
        f_add    = format["add"]
        f_mult   = format["multiply"]
        f_float  = format["floating point"]
        code = {}

        # Loop operands that has to be summed and sort according to map (j,k).
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

            # Exclude all zero valued terms from sum
            value = [v for v in val if not v is None]

            if len(value) > 1:
                # NOTE: Since we no longer call expand_indices, the following
                # is needed to prevent the code from exploding for forms like
                # HyperElasticity
                duplications = {}
                for val in value:
                    if val in duplications:
                        duplications[val] += 1
                        continue
                    duplications[val] = 1

                # Add a product for eacht term that has duplicate code
                expressions = []
                for expr, num_occur in duplications.items():
                    if num_occur > 1:
                        # Pre-multiply expression with number of occurrences
                        expressions.append(f_mult([f_float(num_occur), expr]))
                        continue
                    # Just add expression if there is only one
                    expressions.append(expr)
                ffc_assert(expressions, "Where did the expressions go?")

                if len(expressions) > 1:
                    code[key] = f_group(f_add(expressions))
                    continue
                code[key] = expressions[0]
            else:
                # Check for zero valued sum
                if not value:
                    code[key] = None
                    continue
                code[key] = value[0]

        return code

    def product(self, o, *operands):
        #print("Visiting Product with operands: \n" + "\n".join(map(repr,operands)))

        # Prefetch formats to speed up code generation.
        f_mult = format["multiply"]
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

        #print("\npermute: " + repr(permute))
        #print("\nnot_permute: " + repr(not_permute))
        #print("\npermutations: " + repr(permutations))

        # Create code.
        code ={}
        if permutations:
            for key, val in permutations.items():
                # Sort key in order to create a unique key.
                l = list(key)
                l.sort()

                # Loop products, don't multiply by '1' and if we encounter a None the product is zero.
                # TODO: Need to find a way to remove and J_inv00 terms that might
                # disappear as a consequence of eliminating a zero valued term
                value = []
                zero = False
                for v in val + not_permute:
                    if v is None:
                        ffc_assert(tuple(l) not in code, "This key should not be in the code.")
                        code[tuple(l)] = None
                        zero = True
                        break
                    elif not v:
                        print "v: '%s'" % repr(v)
                        error("should not happen")
                    elif v == "1":
                        pass
                    else:
                        value.append(v)

                if not value:
                    value = ["1"]
                if zero:
                    code[tuple(l)] = None
                else:
                    code[tuple(l)] = f_mult(value)
        else:
            # Loop products, don't multiply by '1' and if we encounter a None the product is zero.
            # TODO: Need to find a way to remove terms from 'used sets' that might
            # disappear as a consequence of eliminating a zero valued term
            value = []
            for v in not_permute:
                if v is None:
                    code[()] = None
                    return code
                elif not v:
                    print "v: '%s'" % repr(v)
                    error("should not happen")
                elif v == "1":
                    pass
                else:
                    value.append(v)
            if value == []:
                value = ["1"]

            code[()] = f_mult(value)

        return code

    def division(self, o, *operands):
        #print("\n\nVisiting Division: " + repr(o) + "with operands: " + "\n".join(map(repr,operands)))

        # Prefetch formats to speed up code generation.
        f_div      = format["div"]
        f_grouping = format["grouping"]

        ffc_assert(len(operands) == 2, \
                   "Expected exactly two operands (numerator and denominator): " + repr(operands))

        # Get the code from the operands.
        numerator_code, denominator_code = operands

        # TODO: Are these safety checks needed? Need to check for None?
        ffc_assert(() in denominator_code and len(denominator_code) == 1, \
                   "Only support function type denominator: " + repr(denominator_code))

        code = {}
        # Get denominator and create new values for the numerator.
        denominator = denominator_code[()]
        ffc_assert(denominator is not None, "Division by zero!")

        for key, val in numerator_code.items():
            # If numerator is None the fraction is also None
            if val is None:
                code[key] = None
            # If denominator is '1', just return numerator
            elif denominator == "1":
                code[key] = val
            # Create fraction and add to code
            else:
                code[key] = f_div(val, f_grouping(denominator))

        return code

    def power(self, o):
        #print("\n\nVisiting Power: " + repr(o))

        # Get base and exponent.
        base, expo = o.operands()

        # Visit base to get base code.
        base_code = self.visit(base)

        # TODO: Are these safety checks needed? Need to check for None?
        ffc_assert(() in base_code and len(base_code) == 1, "Only support function type base: " + repr(base_code))

        # Get the base code.
        val = base_code[()]

        # Handle different exponents
        if isinstance(expo, IntValue):
            return {(): format["power"](val, expo.value())}
        elif isinstance(expo, FloatValue):
            return {(): format["std power"](val, format["floating point"](expo.value()))}
        elif isinstance(expo, (Coefficient, Operator)):
            exp = self.visit(expo)
            return {(): format["std power"](val, exp[()])}
        else:
            error("power does not support this exponent: " + repr(expo))

    def abs(self, o, *operands):
        #print("\n\nVisiting Abs: " + repr(o) + "with operands: " + "\n".join(map(repr,operands)))

        # Prefetch formats to speed up code generation.
        f_abs = format["absolute value"]

        # TODO: Are these safety checks needed? Need to check for None?
        ffc_assert(len(operands) == 1 and () in operands[0] and len(operands[0]) == 1, \
                   "Abs expects one operand of function type: " + repr(operands))

        # Take absolute value of operand.
        return {():f_abs(operands[0][()])}

    # -------------------------------------------------------------------------
    # Condition, Conditional (conditional.py).
    # -------------------------------------------------------------------------
    def not_condition(self, o, *operands):
        # This is a Condition but not a BinaryCondition, and the operand will be another Condition
        # Get condition expression and do safety checks.
        # Might be a bit too strict?
        cond, = operands
        ffc_assert(len(cond) == 1 and cond.keys()[0] == (),\
            "Condition for NotCondition should only be one function: " + repr(cond))
        return {(): format["not"](cond[()])}

    def binary_condition(self, o, *operands):

        # Get LHS and RHS expressions and do safety checks.
        # Might be a bit too strict?
        lhs, rhs = operands
        ffc_assert(len(lhs) == 1 and lhs.keys()[0] == (),\
            "LHS of Condition should only be one function: " + repr(lhs))
        ffc_assert(len(rhs) == 1 and rhs.keys()[0] == (),\
            "RHS of Condition should only be one function: " + repr(rhs))

        # Map names from UFL to cpp.py.
        name_map = {"==":"is equal", "!=":"not equal",\
                    "<":"less than", ">":"greater than",\
                    "<=":"less equal", ">=":"greater equal",\
                    "&&":"and", "||": "or"}

        # Get values and test for None
        l_val = lhs[()]
        r_val = rhs[()]
        if l_val is None:
            l_val = format["float"](0.0)
        if r_val is None:
            r_val = format["float"](0.0)

        return {(): l_val + format[name_map[o._name]] + r_val}

    def conditional(self, o, *operands):

        # Get condition and return values; and do safety check.
        cond, true, false = operands
        ffc_assert(len(cond) == 1 and cond.keys()[0] == (),\
            "Condtion should only be one function: " + repr(cond))
        ffc_assert(len(true) == 1 and true.keys()[0] == (),\
            "True value of Condtional should only be one function: " + repr(true))
        ffc_assert(len(false) == 1 and false.keys()[0] == (),\
            "False value of Condtional should only be one function: " + repr(false))

        # Get values and test for None
        t_val = true[()]
        f_val = false[()]
        if t_val is None:
            t_val = format["float"](0.0)
        if f_val is None:
            f_val = format["float"](0.0)

        # Create expression for conditional
        expr = format["evaluate conditional"](cond[()], t_val, f_val)
        num = len(self.conditionals)
        name = format["conditional"](num)
        if not expr in self.conditionals:
            self.conditionals[expr] = (IP, operation_count(expr, format), num)
        else:
            num = self.conditionals[expr][2]
            name = format["conditional"](num)
        return {():name}

    # -------------------------------------------------------------------------
    # FacetNormal, CellVolume (geometry.py).
    # -------------------------------------------------------------------------
    def facet_normal(self, o,  *operands):
        #print("Visiting FacetNormal:")

        # Get the component
        components = self.component()

        # Safety check.
        ffc_assert(not operands, "Didn't expect any operands for FacetNormal: " + repr(operands))

        # Handle 1D as a special case.
        # FIXME: KBO: This has to change for mD elements in R^n : m < n
        if self.geo_dim == 1:
            # Safety check.
            ffc_assert(len(components) == 0, "FacetNormal in 1D does not expect a component index: " + repr(components))
            normal_component = format["normal component"](self.restriction, "")
            self.trans_set.add(normal_component)
        else:

            # Safety check.
            ffc_assert(len(components) == 1, "FacetNormal expects 1 component index: " + repr(components))

            # We get one component.
            normal_component = format["normal component"](self.restriction, components[0])
            self.trans_set.add(normal_component)

        return {():normal_component}

    def cell_volume(self, o,  *operands):
        # Safety check.
        ffc_assert(not operands, "Didn't expect any operands for CellVolume: " + repr(operands))

        # FIXME: KBO: This has to change for higher order elements
        volume = format["cell volume"](self.restriction)
        self.trans_set.add(volume)

        return {():volume}

    def circumradius(self, o,  *operands):
        # Safety check.
        ffc_assert(not operands, "Didn't expect any operands for Circumradius: " + repr(operands))

        # FIXME: KBO: This has to change for higher order elements
        circumradius = format["circumradius"](self.restriction)
        self.trans_set.add(circumradius)

        return {():circumradius}


    def create_argument(self, ufl_argument, derivatives, component, local_comp,
                  local_offset, ffc_element, transformation, multiindices):
        "Create code for basis functions, and update relevant tables of used basis."

        # Prefetch formats to speed up code generation.
        f_group         = format["grouping"]
        f_add           = format["add"]
        f_mult          = format["multiply"]
        f_transform     = format["transform"]
        f_detJ          = format["det(J)"]
        f_inv           = format["inverse"]

        code = {}
        # Handle affine mappings.
        if transformation == "affine":
            # Loop derivatives and get multi indices.
            for multi in multiindices:
                deriv = [multi.count(i) for i in range(self.geo_dim)]
                if not any(deriv):
                    deriv = []
                # Call function to create mapping and basis name.
                mapping, basis = self._create_mapping_basis(component, deriv, ufl_argument, ffc_element)
                if basis is None:
                    if not mapping in code:
                        code[mapping] = []
                    continue

                # Add transformation if needed.
                if mapping in code:
                    code[mapping].append(self.__apply_transform(basis, derivatives, multi))
                else:
                    code[mapping] = [self.__apply_transform(basis, derivatives, multi)]

        # Handle non-affine mappings.
        else:
            # Loop derivatives and get multi indices.
            for multi in multiindices:
                deriv = [multi.count(i) for i in range(self.geo_dim)]
                if not any(deriv):
                    deriv = []
                for c in range(self.geo_dim):
                    # Call function to create mapping and basis name.
                    mapping, basis = self._create_mapping_basis(c + local_offset, deriv, ufl_argument, ffc_element)
                    if basis is None:
                        if not mapping in code:
                            code[mapping] = []
                        continue

                    # Multiply basis by appropriate transform.
                    if transformation == "covariant piola":
                        dxdX = f_transform("JINV", c, local_comp, self.restriction)
                        self.trans_set.add(dxdX)
                        basis = f_mult([dxdX, basis])
                    elif transformation == "contravariant piola":
                        self.trans_set.add(f_detJ(self.restriction))
                        detJ = f_inv(f_detJ(self.restriction))
                        dXdx = f_transform("J", local_comp, c, self.restriction)
                        self.trans_set.add(dXdx)
                        basis = f_mult([detJ, dXdx, basis])
                    else:
                        error("Transformation is not supported: " + repr(transformation))

                    # Add transformation if needed.
                    if mapping in code:
                        code[mapping].append(self.__apply_transform(basis, derivatives, multi))
                    else:
                        code[mapping] = [self.__apply_transform(basis, derivatives, multi)]

        # Add sums and group if necessary.
        for key, val in code.items():
            if len(val) > 1:
                code[key] = f_group(f_add(val))
            elif val:
                code[key] = val[0]
            else:
                # Return a None (zero) because val == []
                code[key] = None

        return code

    def create_function(self, ufl_function, derivatives, component, local_comp,
                  local_offset, ffc_element, quad_element, transformation, multiindices):
        "Create code for basis functions, and update relevant tables of used basis."

        # Prefetch formats to speed up code generation.
        f_mult          = format["multiply"]
        f_transform     = format["transform"]
        f_detJ          = format["det(J)"]
        f_inv           = format["inverse"]

        code = []
        # Handle affine mappings.
        if transformation == "affine":
            # Loop derivatives and get multi indices.
            for multi in multiindices:
                deriv = [multi.count(i) for i in range(self.geo_dim)]
                if not any(deriv):
                    deriv = []
                # Call other function to create function name.
                function_name = self._create_function_name(component, deriv, quad_element, ufl_function, ffc_element)
                if function_name is None:
                    continue

                # Add transformation if needed.
                code.append(self.__apply_transform(function_name, derivatives, multi))

        # Handle non-affine mappings.
        else:
            # Loop derivatives and get multi indices.
            for multi in multiindices:
                deriv = [multi.count(i) for i in range(self.geo_dim)]
                if not any(deriv):
                    deriv = []
                for c in range(self.geo_dim):
                    function_name = self._create_function_name(c + local_offset, deriv, quad_element, ufl_function, ffc_element)
                    if function_name is None:
                        continue

                    # Multiply basis by appropriate transform.
                    if transformation == "covariant piola":
                        dxdX = f_transform("JINV", c, local_comp, self.restriction)
                        self.trans_set.add(dxdX)
                        function_name = f_mult([dxdX, function_name])
                    elif transformation == "contravariant piola":
                        self.trans_set.add(f_detJ(self.restriction))
                        detJ = f_inv(f_detJ(self.restriction))
                        dXdx = f_transform("J", local_comp, c, self.restriction)
                        self.trans_set.add(dXdx)
                        function_name = f_mult([detJ, dXdx, function_name])
                    else:
                        error("Transformation is not supported: ", repr(transformation))

                    # Add transformation if needed.
                    code.append(self.__apply_transform(function_name, derivatives, multi))

        if not code:
            return None
        elif len(code) > 1:
            code = format["grouping"](format["add"](code))
        else:
            code = code[0]

        return code

    # -------------------------------------------------------------------------
    # Helper functions for Argument and Coefficient
    # -------------------------------------------------------------------------
    def __apply_transform(self, function, derivatives, multi):
        "Apply transformation (from derivatives) to basis or function."
        f_mult          = format["multiply"]
        f_transform     = format["transform"]

        # Add transformation if needed.
        transforms = []
        for i, direction in enumerate(derivatives):
            ref = multi[i]
            t = f_transform("JINV", ref, direction, self.restriction)
            self.trans_set.add(t)
            transforms.append(t)

        # Only multiply by basis if it is present.
        if function:
            prods = transforms + [function]
        else:
            prods = transforms

        return format["multiply"](prods)

    # -------------------------------------------------------------------------
    # Helper functions for transformation of UFL objects in base class
    # -------------------------------------------------------------------------
    def _create_symbol(self, symbol, domain):
        return {():symbol}

    def _create_product(self, symbols):
        return format["multiply"](symbols)

    def _format_scalar_value(self, value):
        #print("format_scalar_value: %d" % value)
        if value is None:
            return {():None}
        # TODO: Handle value < 0 better such that we don't have + -2 in the code.
        return {():format["floating point"](value)}

    def _math_function(self, operands, format_function):
        # TODO: Are these safety checks needed?
        ffc_assert(len(operands) == 1 and () in operands[0] and len(operands[0]) == 1, \
                   "MathFunctions expect one operand of function type: " + repr(operands))
        # Use format function on value of operand.
        new_operand = {}
        operand = operands[0]
        for key, val in operand.items():
            new_operand[key] = format_function(val)
        return new_operand

    # -------------------------------------------------------------------------
    # Helper functions for code_generation()
    # -------------------------------------------------------------------------
    def _count_operations(self, expression):
        return operation_count(expression, format)

    def _create_entry_data(self, val):
        # Multiply value by weight and determinant
        # Create weight and scale factor.
        weight = format["weight"](self.points)
        if self.points > 1:
            weight += format["component"]("", format["integration points"])
        f_scale_factor = format["scale factor"]

        # Update sets of used variables.
        trans_set = set([f_scale_factor])
        trans_set.update(self.trans_set)
        used_points = set([self.points])
        value = format["mul"]([val, weight, f_scale_factor])
        ops = self._count_operations(value)
        used_psi_tables = set([v for k, v in self.psi_tables_map.items()])

        return (value, ops, [trans_set, used_points, used_psi_tables])
