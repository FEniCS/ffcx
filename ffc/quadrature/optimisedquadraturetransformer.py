# -*- coding: utf-8 -*-
"QuadratureTransformer (optimised) for quadrature code generation to translate UFL expressions."

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
# Modified by Anders Logg, 2009

# UFL common.
from ufl.utils.sorting import sorted_by_key
from ufl.measure import custom_integral_types, point_integral_types

# UFL Classes.
from ufl.classes import IntValue
from ufl.classes import FloatValue
from ufl.classes import Coefficient
from ufl.classes import Operator

# FFC modules.
from ffc.log import error, ffc_assert
from ffc.quadrature.cpp import format
from ffc.quadrature.quadraturetransformerbase import QuadratureTransformerBase
from ffc.quadrature.quadratureutils import create_permutations

# Symbolics functions
from ffc.quadrature.symbolics import (create_float, create_symbol,
                                      create_product, create_sum,
                                      create_fraction, BASIS, IP, GEO)


def firstkey(d):
    return next(iter(d))


class QuadratureTransformerOpt(QuadratureTransformerBase):

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
        code = {}
        # Loop operands that has to be summend.
        for op in operands:
            # If entries does already exist we can add the code,
            # otherwise just dump them in the element tensor.
            for key, val in sorted(op.items()):
                if key in code:
                    code[key].append(val)
                else:
                    code[key] = [val]

        # Add sums and group if necessary.
        for key, val in sorted_by_key(code):
            if len(val) > 1:
                code[key] = create_sum(val)
            elif val:
                code[key] = val[0]
            else:
                error("Where did the values go?")
            # If value is zero just ignore it.
            if abs(code[key].val) < format["epsilon"]:
                del code[key]

        return code

    def product(self, o, *operands):

        permute = []
        not_permute = []

        # Sort operands in objects that needs permutation and objects
        # that does not.
        for op in operands:
            # If we get an empty dict, something was zero and so is
            # the product.
            if not op:
                return {}
            if len(op) > 1 or (op and firstkey(op) != ()):
                permute.append(op)
            elif op and firstkey(op) == ():
                not_permute.append(op[()])

        # Create permutations.
        # TODO: After all indices have been expanded I don't think
        # that we'll ever get more than a list of entries and values.
        permutations = create_permutations(permute)

        # Create code.
        code = {}
        if permutations:
            for key, val in permutations.items():
                # Sort key in order to create a unique key.
                l = sorted(key)  # noqa: E741

                # TODO: I think this check can be removed for speed
                # since we just have a list of objects we should never
                # get any conflicts here.
                ffc_assert(tuple(l) not in code,
                           "This key should not be in the code.")

                code[tuple(l)] = create_product(val + not_permute)
        else:
            return {(): create_product(not_permute)}
        return code

    def division(self, o, *operands):

        ffc_assert(len(operands) == 2,
                   "Expected exactly two operands (numerator and denominator): " + repr(operands))

        # Get the code from the operands.
        numerator_code, denominator_code = operands

        # TODO: Are these safety checks needed?
        ffc_assert(() in denominator_code and len(denominator_code) == 1,
                   "Only support function type denominator: " + repr(denominator_code))

        code = {}
        # Get denominator and create new values for the numerator.
        denominator = denominator_code[()]
        for key, val in numerator_code.items():
            code[key] = create_fraction(val, denominator)

        return code

    def power(self, o):
        # Get base and exponent.
        base, expo = o.ufl_operands

        # Visit base to get base code.
        base_code = self.visit(base)

        # TODO: Are these safety checks needed?
        ffc_assert(() in base_code and len(base_code) == 1,
                   "Only support function type base: " + repr(base_code))

        # Get the base code and create power.
        val = base_code[()]

        # Handle different exponents
        if isinstance(expo, IntValue):
            return {(): create_product([val] * expo.value())}
        elif isinstance(expo, FloatValue):
            exp = format["floating point"](expo.value())
            sym = create_symbol(format["std power"](str(val), exp), val.t,
                                val, 1)
            return {(): sym}
        elif isinstance(expo, (Coefficient, Operator)):
            exp = self.visit(expo)[()]
            sym = create_symbol(format["std power"](str(val), exp), val.t,
                                val, 1)
            return {(): sym}
        else:
            error("power does not support this exponent: " + repr(expo))

    def abs(self, o, *operands):

        # TODO: Are these safety checks needed?
        ffc_assert(len(operands) == 1 and () in operands[0] and len(operands[0]) == 1,
                   "Abs expects one operand of function type: " + repr(operands))

        # Take absolute value of operand.
        val = operands[0][()]
        new_val = create_symbol(format["absolute value"](str(val)), val.t,
                                val, 1)
        return {(): new_val}

    def min_value(self, o, *operands):
        # Take minimum value of operands.
        val0 = operands[0][()]
        val1 = operands[1][()]
        t = min(val0.t, val1.t)
        # FIXME: I don't know how to implement this the optimized
        # way. Is this right?
        new_val = create_symbol(format["min value"](str(val0), str(val1)), t)
        return {(): new_val}

    def max_value(self, o, *operands):
        # Take maximum value of operands.
        val0 = operands[0][()]
        val1 = operands[1][()]
        t = min(val0.t, val1.t)
        # FIXME: I don't know how to implement this the optimized
        # way. Is this right?
        new_val = create_symbol(format["max value"](str(val0), str(val1)), t)
        return {(): new_val}

    # -------------------------------------------------------------------------
    # Condition, Conditional (conditional.py).
    # -------------------------------------------------------------------------
    def not_condition(self, o, *operands):
        # This is a Condition but not a BinaryCondition, and the
        # operand will be another Condition
        # Get condition expression and do safety checks.
        # Might be a bit too strict?
        c, = operands
        ffc_assert(len(c) == 1 and firstkey(c) == (),
                   "Condition for NotCondition should only be one function: " + repr(c))
        sym = create_symbol(format["not"](str(c[()])), c[()].t, base_op=c[()].ops() + 1)
        return {(): sym}

    def binary_condition(self, o, *operands):

        # Get LHS and RHS expressions and do safety checks.  Might be
        # a bit too strict?
        lhs, rhs = operands
        ffc_assert(len(lhs) == 1 and firstkey(lhs) == (),
                   "LHS of Condtion should only be one function: " + repr(lhs))
        ffc_assert(len(rhs) == 1 and firstkey(rhs) == (),
                   "RHS of Condtion should only be one function: " + repr(rhs))

        # Map names from UFL to cpp.py.
        name_map = {"==": "is equal", "!=": "not equal",
                    "<": "less than", ">": "greater than",
                    "<=": "less equal", ">=": "greater equal",
                    "&&": "and", "||": "or"}

        # Get the minimum type
        t = min(lhs[()].t, rhs[()].t)
        ops = lhs[()].ops() + rhs[()].ops() + 1
        cond = str(lhs[()]) + format[name_map[o._name]] + str(rhs[()])
        sym = create_symbol(format["grouping"](cond), t, base_op=ops)
        return {(): sym}

    def conditional(self, o, *operands):
        # Get condition and return values; and do safety check.
        cond, true, false = operands
        ffc_assert(len(cond) == 1 and firstkey(cond) == (),
                   "Condtion should only be one function: " + repr(cond))
        ffc_assert(len(true) == 1 and firstkey(true) == (),
                   "True value of Condtional should only be one function: " + repr(true))
        ffc_assert(len(false) == 1 and firstkey(false) == (),
                   "False value of Condtional should only be one function: " + repr(false))

        # Get values and test for None
        t_val = true[()]
        f_val = false[()]

        # Get the minimum type and number of operations
        # TODO: conditionals are currently always located inside the
        # ip loop, therefore the type has to be at least IP (fix bug
        # #1082048). This can be optimised.
        t = min([cond[()].t, t_val.t, f_val.t, IP])
        ops = sum([cond[()].ops(), t_val.ops(), f_val.ops()])

        # Create expression for conditional
        # TODO: Handle this differently to expose the variables which
        # are used to create the expressions.
        expr = create_symbol(format["evaluate conditional"](cond[()], t_val,
                                                            f_val), t)
        num = len(self.conditionals)
        name = create_symbol(format["conditional"](num), t)
        if expr not in self.conditionals:
            self.conditionals[expr] = (t, ops, num)
        else:
            num = self.conditionals[expr][2]
            name = create_symbol(format["conditional"](num), t)
        return {(): name}

    # -------------------------------------------------------------------------
    # FacetNormal, CellVolume, Circumradius, FacetArea (geometry.py).
    # -------------------------------------------------------------------------
    def cell_coordinate(self, o):  # FIXME
        error("This object should be implemented by the child class.")

    def facet_coordinate(self, o):  # FIXME
        error("This object should be implemented by the child class.")

    def cell_origin(self, o):  # FIXME
        error("This object should be implemented by the child class.")

    def facet_origin(self, o):  # FIXME
        error("This object should be implemented by the child class.")

    def cell_facet_origin(self, o):  # FIXME
        error("This object should be implemented by the child class.")

    def jacobian(self, o):  # FIXME
        error("This object should be implemented by the child class.")

    def jacobian_determinant(self, o):  # FIXME
        error("This object should be implemented by the child class.")

    def jacobian_inverse(self, o):  # FIXME
        error("This object should be implemented by the child class.")

    def facet_jacobian(self, o):  # FIXME
        error("This object should be implemented by the child class.")

    def facet_jacobian_determinant(self, o):  # FIXME
        error("This object should be implemented by the child class.")

    def facet_jacobian_inverse(self, o):  # FIXME
        error("This object should be implemented by the child class.")

    def cell_facet_jacobian(self, o):  # FIXME
        error("This object should be implemented by the child class.")

    def cell_facet_jacobian_determinant(self, o):  # FIXME
        error("This object should be implemented by the child class.")

    def cell_facet_jacobian_inverse(self, o):  # FIXME
        error("This object should be implemented by the child class.")

    def facet_normal(self, o):
        components = self.component()

        # Safety check.
        ffc_assert(len(components) == 1,
                   "FacetNormal expects 1 component index: " + repr(components))

        # Handle 1D as a special case.
        # FIXME: KBO: This has to change for mD elements in R^n : m <
        # n
        if self.gdim == 1:  # FIXME: MSA UFL uses shape (1,) now, can we remove the special case here then?
            normal_component = format["normal component"](self.restriction, "")
        else:
            normal_component = format["normal component"](self.restriction,
                                                          components[0])
        self.trans_set.add(normal_component)

        return {(): create_symbol(normal_component, GEO)}

    def cell_normal(self, o):  # FIXME
        error("This object should be implemented by the child class.")

    def cell_volume(self, o):
        # FIXME: KBO: This has to change for higher order elements
        # detJ = format["det(J)"](self.restriction)
        # volume = format["absolute value"](detJ)
        # self.trans_set.add(detJ)

        volume = format["cell volume"](self.restriction)
        self.trans_set.add(volume)

        return {(): create_symbol(volume, GEO)}

    def circumradius(self, o):
        # FIXME: KBO: This has to change for higher order elements
        circumradius = format["circumradius"](self.restriction)
        self.trans_set.add(circumradius)

        return {(): create_symbol(circumradius, GEO)}

    def facet_area(self, o):
        # FIXME: KBO: This has to change for higher order elements
        # NOTE: Omitting restriction because the area of a facet is
        # the same on both sides.
        # FIXME: Since we use the scale factor, facet area has no
        # meaning for cell integrals. (Need check in FFC or UFL).
        area = format["facet area"]
        self.trans_set.add(area)

        return {(): create_symbol(area, GEO)}

    def min_facet_edge_length(self, o):
        # FIXME: this has no meaning for cell integrals. (Need check
        # in FFC or UFL).

        tdim = self.tdim
        if tdim < 3:
            return self.facet_area(o)

        edgelen = format["min facet edge length"](self.restriction)
        self.trans_set.add(edgelen)

        return {(): create_symbol(edgelen, GEO)}

    def max_facet_edge_length(self, o):
        # FIXME: this has no meaning for cell integrals. (Need check
        # in FFC or UFL).

        tdim = self.tdim
        if tdim < 3:
            return self.facet_area(o)

        edgelen = format["max facet edge length"](self.restriction)
        self.trans_set.add(edgelen)

        return {(): create_symbol(edgelen, GEO)}

    def cell_orientation(self, o):  # FIXME
        error("This object should be implemented by the child class.")

    def quadrature_weight(self, o):  # FIXME
        error("This object should be implemented by the child class.")

    def create_argument(self, ufl_argument, derivatives, component, local_comp,
                        local_offset, ffc_element, transformation,
                        multiindices, tdim, gdim, avg):
        "Create code for basis functions, and update relevant tables of used basis."

        # Prefetch formats to speed up code generation.
        f_transform = format["transform"]
        f_detJ = format["det(J)"]

        # Reset code
        code = {}

        # Affine mapping
        if transformation == "affine":
            # Loop derivatives and get multi indices.
            for multi in multiindices:
                deriv = [multi.count(i) for i in range(tdim)]
                if not any(deriv):
                    deriv = []

                # Create mapping and basis name.
                mapping, basis = self._create_mapping_basis(component, deriv,
                                                            avg, ufl_argument,
                                                            ffc_element)
                if mapping not in code:
                    code[mapping] = []

                if basis is not None:
                    # Add transformation if needed.
                    code[mapping].append(self.__apply_transform(basis,
                                                                derivatives,
                                                                multi, tdim,
                                                                gdim))

        # Handle non-affine mappings.
        else:
            ffc_assert(avg is None,
                       "Taking average is not supported for non-affine mappings.")

            # Loop derivatives and get multi indices.
            for multi in multiindices:
                deriv = [multi.count(i) for i in range(tdim)]
                if not any(deriv):
                    deriv = []

                if transformation in ["covariant piola",
                                      "contravariant piola"]:
                    for c in range(tdim):
                        # Create mapping and basis name.
                        mapping, basis = self._create_mapping_basis(c + local_offset, deriv, avg, ufl_argument, ffc_element)
                        if mapping not in code:
                            code[mapping] = []

                        if basis is not None:
                            # Multiply basis by appropriate transform.
                            if transformation == "covariant piola":
                                dxdX = create_symbol(f_transform("JINV", c,
                                                                 local_comp, tdim,
                                                                 gdim,
                                                                 self.restriction),
                                                     GEO)
                                basis = create_product([dxdX, basis])
                            elif transformation == "contravariant piola":
                                detJ = create_fraction(create_float(1),
                                                       create_symbol(f_detJ(self.restriction), GEO))
                                dXdx = create_symbol(f_transform("J", local_comp,
                                                                 c, gdim, tdim,
                                                                 self.restriction),
                                                     GEO)
                                basis = create_product([detJ, dXdx, basis])
                        # Add transformation if needed.
                        code[mapping].append(self.__apply_transform(basis,
                                                                    derivatives,
                                                                    multi, tdim,
                                                                    gdim))
                elif transformation == "double covariant piola":
                    # g_ij = (Jinv)_ki G_kl (Jinv)lj
                    i = local_comp // tdim
                    j = local_comp % tdim
                    for k in range(tdim):
                        for l in range(tdim):
                            # Create mapping and basis name.
                            mapping, basis = self._create_mapping_basis(
                                k * tdim + l + local_offset,
                                deriv, avg, ufl_argument, ffc_element)
                            if mapping not in code:
                                code[mapping] = []
                            if basis is not None:
                                J1 = create_symbol(
                                    f_transform("JINV", k, i, tdim, gdim,
                                                self.restriction), GEO)
                                J2 = create_symbol(
                                    f_transform("JINV", l, j, tdim, gdim,
                                                 self.restriction), GEO)
                                basis = create_product([J1, basis, J2])
                                # Add transformation if needed.
                                code[mapping].append(
                                    self.__apply_transform(
                                        basis, derivatives, multi,
                                        tdim, gdim))
                elif transformation == "double contravariant piola":
                    # g_ij = (detJ)^(-2) J_ik G_kl J_jl
                    i = local_comp // tdim
                    j = local_comp % tdim
                    for k in range(tdim):
                        for l in range(tdim):
                            # Create mapping and basis name.
                            mapping, basis = self._create_mapping_basis(
                                k * tdim + l + local_offset,
                                deriv, avg, ufl_argument, ffc_element)
                            if mapping not in code:
                                code[mapping] = []
                            if basis is not None:
                                J1 = create_symbol(
                                    f_transform("J", i, k, gdim, tdim,
                                                self.restriction), GEO)
                                J2 = create_symbol(
                                    f_transform("J", j, l, gdim, tdim,
                                                self.restriction), GEO)
                                invdetJ = create_fraction(
                                    create_float(1),
                                    create_symbol(f_detJ(self.restriction),
                                                  GEO))
                                basis = create_product([invdetJ, invdetJ, J1,
                                                        basis, J2])
                                # Add transformation if needed.
                                code[mapping].append(
                                    self.__apply_transform(
                                        basis, derivatives, multi,
                                        tdim, gdim))
                else:
                    error("Transformation is not supported: " + repr(transformation))


        # Add sums and group if necessary.
        for key, val in list(code.items()):
            if len(val) > 1:
                code[key] = create_sum(val)
            else:
                code[key] = val[0]

        return code

    def create_function(self, ufl_function, derivatives, component, local_comp,
                        local_offset, ffc_element, is_quad_element,
                        transformation, multiindices, tdim, gdim, avg):
        "Create code for basis functions, and update relevant tables of used basis."
        ffc_assert(ufl_function in self._function_replace_values,
                   "Expecting ufl_function to have been mapped prior to this call.")

        # Prefetch formats to speed up code generation.
        f_transform = format["transform"]
        f_detJ = format["det(J)"]

        # Reset code
        code = []

        # Handle affine mappings.
        if transformation == "affine":
            # Loop derivatives and get multi indices.
            for multi in multiindices:
                deriv = [multi.count(i) for i in range(tdim)]
                if not any(deriv):
                    deriv = []

                # Create function name.
                function_name = self._create_function_name(component, deriv,
                                                           avg, is_quad_element,
                                                           ufl_function,
                                                           ffc_element)
                if function_name:
                    # Add transformation if needed.
                    code.append(self.__apply_transform(function_name,
                                                       derivatives, multi, tdim, gdim))

        # Handle non-affine mappings.
        else:
            ffc_assert(avg is None,
                       "Taking average is not supported for non-affine mappings.")

            # Loop derivatives and get multi indices.
            for multi in multiindices:
                deriv = [multi.count(i) for i in range(tdim)]
                if not any(deriv):
                    deriv = []

                if transformation in ["covariant piola",
                                      "contravariant piola"]:
                    for c in range(tdim):
                        function_name = self._create_function_name(c + local_offset, deriv, avg, is_quad_element, ufl_function, ffc_element)
                        if function_name:
                            # Multiply basis by appropriate transform.
                            if transformation == "covariant piola":
                                dxdX = create_symbol(f_transform("JINV", c,
                                                                 local_comp, tdim,
                                                                 gdim,
                                                                 self.restriction),
                                                     GEO)
                                function_name = create_product([dxdX, function_name])
                            elif transformation == "contravariant piola":
                                detJ = create_fraction(create_float(1),
                                                       create_symbol(f_detJ(self.restriction),
                                                                     GEO))
                                dXdx = create_symbol(f_transform("J", local_comp,
                                                                 c, gdim, tdim,
                                                                 self.restriction),
                                                     GEO)
                                function_name = create_product([detJ, dXdx,
                                                                function_name])
                            # Add transformation if needed.
                            code.append(self.__apply_transform(function_name,
                                                               derivatives, multi,
                                                               tdim, gdim))
                elif transformation == "double covariant piola":
                    # g_ij = (Jinv)_ki G_kl (Jinv)lj
                    i = local_comp // tdim
                    j = local_comp % tdim
                    for k in range(tdim):
                        for l in range(tdim):
                            # Create mapping and basis name.
                            function_name = self._create_function_name(
                                k * tdim + l + local_offset, deriv, avg,
                                is_quad_element, ufl_function, ffc_element)
                            J1 = create_symbol(
                                f_transform("JINV", k, i, tdim, gdim,
                                            self.restriction), GEO)
                            J2 = create_symbol(
                                f_transform("JINV", l, j, tdim, gdim,
                                            self.restriction), GEO)
                            function_name = create_product([J1, function_name,
                                                            J2])
                            # Add transformation if needed.
                            code.append(self.__apply_transform(
                                function_name, derivatives, multi, tdim, gdim))
                elif transformation == "double contravariant piola":
                    # g_ij = (detJ)^(-2) J_ik G_kl J_jl
                    i = local_comp // tdim
                    j = local_comp % tdim
                    for k in range(tdim):
                        for l in range(tdim):
                            # Create mapping and basis name.
                            function_name = self._create_function_name(
                                k * tdim + l + local_offset,
                                deriv, avg, is_quad_element,
                                ufl_function, ffc_element)
                            J1 = create_symbol(
                                f_transform("J", i, k, tdim, gdim,
                                            self.restriction), GEO)
                            J2 = create_symbol(
                                f_transform("J", j, l, tdim, gdim,
                                            self.restriction), GEO)
                            invdetJ = create_fraction(
                                create_float(1),
                                create_symbol(f_detJ(self.restriction), GEO))
                            function_name = create_product([invdetJ, invdetJ,
                                                            J1, function_name,
                                                            J2])
                            # Add transformation if needed.
                            code.append(self.__apply_transform(function_name,
                                                               derivatives,
                                                               multi, tdim,
                                                               gdim))
                else:
                    error("Transformation is not supported: ",
                          repr(transformation))

        if not code:
            return create_float(0.0)
        elif len(code) > 1:
            code = create_sum(code)
        else:
            code = code[0]

        return code

    # -------------------------------------------------------------------------
    # Helper functions for Argument and Coefficient
    # -------------------------------------------------------------------------
    def __apply_transform(self, function, derivatives, multi, tdim, gdim):
        "Apply transformation (from derivatives) to basis or function."
        f_transform = format["transform"]

        # Add transformation if needed.
        transforms = []
        if self.integral_type not in custom_integral_types:
            for i, direction in enumerate(derivatives):
                ref = multi[i]
                t = f_transform("JINV", ref, direction, tdim, gdim,
                                self.restriction)
                transforms.append(create_symbol(t, GEO))

        transforms.append(function)
        return create_product(transforms)

    # -------------------------------------------------------------------------
    # Helper functions for transformation of UFL objects in base class
    # -------------------------------------------------------------------------
    def _create_symbol(self, symbol, domain):
        return {(): create_symbol(symbol, domain)}

    def _create_product(self, symbols):
        return create_product(symbols)

    def _format_scalar_value(self, value):
        # print("format_scalar_value: %d" % value)
        if value is None:
            return {(): create_float(0.0)}
        return {(): create_float(value)}

    def _math_function(self, operands, format_function):
        # TODO: Are these safety checks needed?
        ffc_assert(len(operands) == 1 and () in operands[0] and len(operands[0]) == 1,
                   "MathFunctions expect one operand of function type: " + repr(operands))
        # Use format function on value of operand.
        operand = operands[0]
        for key, val in list(operand.items()):
            new_val = create_symbol(format_function(str(val)), val.t, val, 1)
            operand[key] = new_val
        # raise Exception("pause")
        return operand

    def _bessel_function(self, operands, format_function):
        # TODO: Are these safety checks needed?
        # TODO: work on reference instead of copies? (like
        # math_function)
        ffc_assert(len(operands) == 2,
                   "BesselFunctions expect two operands of function type: " + repr(operands))
        nu, x = operands
        ffc_assert(len(nu) == 1 and () in nu,
                   "Expecting one operand of function type as first argument to BesselFunction : " + repr(nu))
        ffc_assert(len(x) == 1 and () in x,
                   "Expecting one operand of function type as second argument to BesselFunction : " + repr(x))
        nu = nu[()]
        x = x[()]
        if nu is None:
            nu = format["floating point"](0.0)
        if x is None:
            x = format["floating point"](0.0)

        sym = create_symbol(format_function(x, nu), x.t, x, 1)
        return {(): sym}

    # -------------------------------------------------------------------------
    # Helper functions for code_generation()
    # -------------------------------------------------------------------------
    def _count_operations(self, expression):
        return expression.ops()

    def _create_entry_data(self, val, integral_type):
        # Multiply value by weight and determinant
        ACCESS = GEO
        weight = format["weight"](self.points)
        if self.points > 1:
            weight += format["component"]("", format["integration points"])
            ACCESS = IP
        weight = self._create_symbol(weight, ACCESS)[()]

        # Create value.
        if integral_type in (point_integral_types + custom_integral_types):
            trans_set = set()
            value = create_product([val, weight])
        else:
            f_scale_factor = format["scale factor"]
            trans_set = set([f_scale_factor])
            value = create_product([val, weight,
                                    create_symbol(f_scale_factor, GEO)])

        # Update sets of used variables (if they will not be used
        # because of optimisations later, they will be reset).
        trans_set.update([str(x) for x in value.get_unique_vars(GEO)])
        used_points = set([self.points])
        ops = self._count_operations(value)
        used_psi_tables = set([self.psi_tables_map[b]
                               for b in value.get_unique_vars(BASIS)])

        return (value, ops, [trans_set, used_points, used_psi_tables])
