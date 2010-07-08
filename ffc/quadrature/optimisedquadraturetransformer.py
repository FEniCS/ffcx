"QuadratureTransformer (optimised) for quadrature code generation to translate UFL expressions."

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@gmail.com)"
__date__ = "2009-03-18"
__copyright__ = "Copyright (C) 2009-2010 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Anders Logg, 2009
# Last changed: 2010-03-11

# Python modules.
from numpy import shape

# UFL common.
from ufl.common import product

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
from ffc.quadrature.quadraturetransformerbase import QuadratureTransformerBase
from ffc.quadrature.quadratureutils import create_permutations

# Symbolics functions
#from symbolics import set_format
from ffc.quadrature.symbolics import create_float, create_symbol, create_product,\
                                     create_sum, create_fraction, BASIS, IP, GEO, CONST

class QuadratureTransformerOpt(QuadratureTransformerBase):
    "Transform UFL representation to quadrature code."

    def __init__(self, *args):

        # Initialise base class.
        QuadratureTransformerBase.__init__(self, *args)
#        set_format(format)

    # -------------------------------------------------------------------------
    # Start handling UFL classes.
    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------
    # AlgebraOperators (algebra.py).
    # -------------------------------------------------------------------------
    def sum(self, o, *operands):
        #print("Visiting Sum: " + repr(o) + "\noperands: " + "\n".join(map(repr, operands)))

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
                code[key] = create_sum(val)
            elif val:
                code[key] = val[0]
            else:
                error("Where did the values go?")

        return code

    def product(self, o, *operands):
        #print("\n\nVisiting Product:\n" + str(tree_format(o)))

        permute = []
        not_permute = []

        # Sort operands in objects that needs permutation and objects that does not.
        for op in operands:
            if len(op) > 1 or (op and op.keys()[0] != ()):
                permute.append(op)
            elif op:
                not_permute.append(op[()])

        # Create permutations.
        # TODO: After all indices have been expanded I don't think that we'll
        # ever get more than a list of entries and values.
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

                # TODO: I think this check can be removed for speed since we
                # just have a list of objects we should never get any conflicts here.
                ffc_assert(tuple(l) not in code, "This key should not be in the code.")

                code[tuple(l)] = create_product(val + not_permute)
        else:
            return {():create_product(not_permute)}
        return code

    def division(self, o, *operands):
        #print("\n\nVisiting Division: " + repr(o) + "with operands: " + "\n".join(map(repr,operands)))

        ffc_assert(len(operands) == 2, "Expected exactly two operands (numerator and denominator): " + repr(operands))

        # Get the code from the operands.
        numerator_code, denominator_code = operands

        # TODO: Are these safety checks needed?
        ffc_assert(() in denominator_code and len(denominator_code) == 1, \
                   "Only support function type denominator: " + repr(denominator_code))

        code = {}
        # Get denominator and create new values for the numerator.
        denominator = denominator_code[()]
        for key, val in numerator_code.items():
            code[key] = create_fraction(val, denominator)

        return code

    def power(self, o):
        #print("\n\nVisiting Power: " + repr(o))

        # Get base and exponent.
        base, expo = o.operands()

        # Visit base to get base code.
        base_code = self.visit(base)

        # TODO: Are these safety checks needed?
        ffc_assert(() in base_code and len(base_code) == 1, "Only support function type base: " + repr(base_code))

        # Get the base code and create power.
        val = base_code[()]

        # Handle different exponents
        if isinstance(expo, IntValue):
            return {(): create_product([val]*expo.value())}
        elif isinstance(expo, FloatValue):
            exp = format["floating point"](expo.value())
#            sym = create_symbol(format["std power"](str(val), exp), val.t)
#            sym.base_expr = val
#            sym.base_op = 1 # Add one operation for the pow() function.
            sym = create_symbol(format["std power"], val.t, val, 1)
            sym.exp = exp
            return {(): sym}
        elif isinstance(expo, (Coefficient, Operator)):
            exp = self.visit(expo)
#            sym = create_symbol(format["std power"](str(val), exp[()]), val.t)
#            sym.base_expr = val
#            sym.base_op = 1 # Add one operation for the pow() function.
            sym = create_symbol(format["std power"], val.t, val, 1)
            sym.exp = exp[()]
            return {(): sym}
        else:
            error("power does not support this exponent: " + repr(expo))

    def abs(self, o, *operands):
        #print("\n\nVisiting Abs: " + repr(o) + "with operands: " + "\n".join(map(repr,operands)))

        # TODO: Are these safety checks needed?
        ffc_assert(len(operands) == 1 and () in operands[0] and len(operands[0]) == 1, \
                   "Abs expects one operand of function type: " + repr(operands))

        # Take absolute value of operand.
        val = operands[0][()]
#        new_val = create_symbol(format["absolute value"](str(val)), val.t)
#        new_val.base_expr = val
#        new_val.base_op = 1 # Add one operation for taking the absolute value.
        new_val = create_symbol(format["absolute value"], val.t, val, 1)
        return {():new_val}

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

        return {(): create_symbol(normal_component, GEO)}

    def cell_volume(self, o,  *operands):
        # Safety check.
        ffc_assert(not operands, "Didn't expect any operands for FacetNormal: " + repr(operands))

        # FIXME: KBO: This has to change for higher order elements
#        detJ = format["det(J)"](self.restriction)
#        volume = format["absolute value"](detJ)
#        self.trans_set.add(detJ)

        volume = format["cell volume"](self.restriction)
        self.trans_set.add(volume)

        return {():create_symbol(volume, GEO)}

    def circumradius(self, o,  *operands):
        # Safety check.
        ffc_assert(not operands, "Didn't expect any operands for Circumradius: " + repr(operands))

        # FIXME: KBO: This has to change for higher order elements
        circumradius = format["circumradius"](self.restriction)
        self.trans_set.add(circumradius)

        return {():create_symbol(circumradius, GEO)}

    def create_argument(self, ufl_argument, derivatives, component, local_comp,
                  local_offset, ffc_element, transformation, multiindices):
        "Create code for basis functions, and update relevant tables of used basis."

        # Prefetch formats to speed up code generation.
        f_transform     = format["transform"]
        f_detJ          = format["det(J)"]

        code = {}

        # Affince mapping
        if transformation == "affine":
            # Loop derivatives and get multi indices.
            for multi in multiindices:
                deriv = [multi.count(i) for i in range(self.geo_dim)]
                if not any(deriv):
                    deriv = []

                # Call function to create mapping and basis name.
                mapping, basis = self._create_mapping_basis(component, deriv, ufl_argument, ffc_element)

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

                    # Multiply basis by appropriate transform.
                    if transformation == "covariant piola":
                        dxdX = create_symbol(f_transform("JINV", c, local_comp, self.restriction), GEO)
                        basis = create_product([dxdX, basis])
                    elif transformation == "contravariant piola":
                        detJ = create_fraction(create_float(1), create_symbol(f_detJ(self.restriction), GEO))
                        dXdx = create_symbol(f_transform("J", local_comp, c, self.restriction), GEO)
                        basis = create_product([detJ, dXdx, basis])
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
                code[key] = create_sum(val)
            else:
                code[key] = val[0]

        return code

    def create_function(self, ufl_function, derivatives, component, local_comp,
                  local_offset, ffc_element, quad_element, transformation, multiindices):
        "Create code for basis functions, and update relevant tables of used basis."

        # Prefetch formats to speed up code generation.
        f_transform     = format["transform"]
        f_detJ          = format["det(J)"]

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
                if not function_name:
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

                    # Multiply basis by appropriate transform.
                    if transformation == "covariant piola":
                        dxdX = create_symbol(f_transform("JINV", c, local_comp, self.restriction), GEO)
                        function_name = create_product([dxdX, function_name])
                    elif transformation == "contravariant piola":
                        detJ = create_fraction(create_float(1), create_symbol(f_detJ(self.restriction), GEO))
                        dXdx = create_symbol(f_transform("J", local_comp, c, self.restriction), GEO)
                        function_name = create_product([detJ, dXdx, function_name])
                    else:
                        error("Transformation is not supported: ", repr(transformation))

                    # Add transformation if needed.
                    code.append(self.__apply_transform(function_name, derivatives, multi))
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
    def __apply_transform(self, function, derivatives, multi):
        "Apply transformation (from derivatives) to basis or function."
        f_transform     = format["transform"]

        # Add transformation if needed.
        transforms = []
        for i, direction in enumerate(derivatives):
            ref = multi[i]
            t = f_transform("JINV", ref, direction, self.restriction)
            transforms.append(create_symbol(t, GEO))
        transforms.append(function)
        return create_product(transforms)

    # -------------------------------------------------------------------------
    # Helper functions for transformation of UFL objects in base class
    # -------------------------------------------------------------------------
    def _create_symbol(self, symbol, domain):
        return {():create_symbol(symbol, domain)}

    def _create_product(self, symbols):
        return create_product(symbols)

    def _format_scalar_value(self, value):
        #print("format_scalar_value: %d" % value)
        if value is None:
            return {():create_float(0.0)}
        return {():create_float(value)}

    def _math_function(self, operands, format_function):
        #print("Calling _math_function() of optimisedquadraturetransformer.")
        # TODO: Are these safety checks needed?
        ffc_assert(len(operands) == 1 and () in operands[0] and len(operands[0]) == 1, \
                   "MathFunctions expect one operand of function type: " + repr(operands))
        # Use format function on value of operand.
        operand = operands[0]
        for key, val in operand.items():
#            new_val = create_symbol(format_function(str(val)), val.t)
#            new_val.base_expr = val
#            new_val.base_op = 1 # Add one operation for the math function.
            new_val = create_symbol(format_function, val.t, val, 1)
            operand[key] = new_val
        return operand

    # -------------------------------------------------------------------------
    # Helper functions for code_generation()
    # -------------------------------------------------------------------------
    def _count_operations(self, expression):
        return expression.ops()

    def _create_entry_data(self, val):
#        zero = False
        # Multiply value by weight and determinant
        ACCESS = GEO
        weight = format["weight"](self.points)
        if self.points > 1:
            weight += format["component"]("", format["integration points"])
            ACCESS = IP
        weight = self._create_symbol(weight, ACCESS)[()]
        f_scale_factor = format["scale factor"]

        # Create value.
        value = create_product([val, weight, create_symbol(f_scale_factor, GEO)])

        # Update sets of used variables (if they will not be used because of
        # optimisations later, they will be reset).
        trans_set = set([f_scale_factor])
        trans_set.update(map(lambda x: str(x), value.get_unique_vars(GEO)))
        used_points = set([self.points])
        ops = self._count_operations(value)
        used_psi_tables = set([self.psi_tables_map[b] for b in value.get_unique_vars(BASIS)])

        return (value, ops, [trans_set, used_points, used_psi_tables])

