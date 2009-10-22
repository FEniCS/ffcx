"QuadratureTransformer for quadrature code generation to translate UFL expressions."

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2009-02-09 -- 2009-10-19"
__copyright__ = "Copyright (C) 2009 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Peter Brune 2009

# Python modules.
from numpy import shape

# UFL common.
from ufl.common import product, StackDict, Stack

# UFL Classes.
from ufl.classes import FixedIndex
from ufl.classes import IntValue
from ufl.classes import FloatValue
from ufl.classes import Function

# UFL Algorithms.
from ufl.algorithms.printing import tree_format

# FFC common modules.
from ffc.common.log import info, debug, error

# FFC fem modules.
from ffc.fem.finiteelement import AFFINE, CONTRAVARIANT_PIOLA, COVARIANT_PIOLA

# Utility and optimisation functions for quadraturegenerator.
from quadraturetransformerbase import QuadratureTransformerBase
from quadraturegenerator_utils import generate_psi_name
from quadraturegenerator_utils import create_permutations
from reduce_operations import operation_count

class QuadratureTransformer(QuadratureTransformerBase):
    "Transform UFL representation to quadrature code."

    def __init__(self, form_representation, domain_type, optimise_options, format):

        QuadratureTransformerBase.__init__(self, form_representation, domain_type, optimise_options, format)

    # -------------------------------------------------------------------------
    # Start handling UFL classes.
    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------
    # AlgebraOperators (algebra.py).
    # -------------------------------------------------------------------------
    def sum(self, o, *operands):
        #print("Visiting Sum: " + "\noperands: \n" + "\n".join(map(str, operands)))

        # Prefetch formats to speed up code generation.
        format_group  = self.format["grouping"]
        format_add    = self.format["add"]
        format_mult   = self.format["multiply"]
        format_float  = self.format["floating point"]
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
                        expressions.append(format_mult([format_float(num_occur), expr]))
                        continue
                    # Just add expression if there is only one
                    expressions.append(expr)

                if not expressions:
                    error("Where did the expressions go?")
                if len(expressions) > 1:
                    code[key] = format_group(format_add(expressions))
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
        #print("Visiting Product with operands: \n" + "\n".join(map(str,operands)))

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

        #print("\npermute: " + str(permute))
        #print("\nnot_permute: " + str(not_permute))
        #print("\npermutations: " + str(permutations))

        # Create code.
        code ={}
        if permutations:
            for key, val in permutations.items():
                # Sort key in order to create a unique key.
                l = list(key)
                l.sort()

#                value = [v for v in val + not_permute if v]
#                value = [v for v in val + not_permute]
#                value = [v for v in val + not_permute if v and v != "1"]

                # Loop products, don't multiply by '1' and if we encounter a None the product is zero.
                # TODO: Need to find a way to remove and J_inv00 terms that might
                # disappear as a consequence of eliminating a zero valued term
                value = []
                zero = False
                for v in val + not_permute:
                    if v is None:
                        #print "v is None"
                        if tuple(l) in code:
                            error("This key should not be in the code.")
                        code[tuple(l)] = None
                        zero = True
                        break
                    elif not v:
                        print "v: '%s'" % str(v)
                        raise RuntimeError("should not happen")
                    elif v == "1":
                        #print "v == 1:"
                        pass
                    else:
                        value.append(v)

                if not value:
                    #print "No values left: '%s'" % format_mult(value)
                    value = ["1"]
                if zero:
                    code[tuple(l)] = None
                else:
                    code[tuple(l)] = format_mult(value)
        else:
#            value = [v for v in not_permute if v]
#            value = [v for v in not_permute if v and v != "1"]

            # Loop products, don't multiply by '1' and if we encounter a None the product is zero.
            # TODO: Need to find a way to remove terms from 'used sets' that might
            # disappear as a consequence of eliminating a zero valued term
            value = []
            for v in not_permute:
                if v is None:
                    #print "v is None"
                    code[()] = None
                    return code
                elif not v:
                    print "v: '%s'" % str(v)
                    raise RuntimeError("should not happen")
                elif v == "1":
                    #print "v == 1: '%s'" % str(v)
                    pass
                else:
                    value.append(v)

            if value == []:
                #print "No values left: '%s'" % format_mult(value)
                value = ["1"]

            code[()] = format_mult(value)

        return code

    def division(self, o, *operands):
        #print("\n\nVisiting Division: " + o.__repr__() + "with operands: " + "\n".join(map(str,operands)))

        # Prefetch formats to speed up code generation.
        format_div      = self.format["division"]
        format_grouping = self.format["grouping"]

        if len(operands) != 2:
            error("Expected exactly two operands (numerator and denominator): " + operands.__repr__())

        # Get the code from the operands.
        numerator_code, denominator_code = operands
        #print("\nnumerator: " + str(numerator_code))
        #print("\ndenominator: " + str(denominator_code))

        # TODO: Are these safety checks needed? Need to check for None?
        if not () in denominator_code and len(denominator_code) != 1:
            error("Only support function type denominator: " + str(denominator_code))

        code = {}
        # Get denominator and create new values for the numerator.
        denominator = denominator_code[()]
        for key, val in numerator_code.items():
            code[key] = val + format_div + format_grouping(denominator)

        return code

    def power(self, o):
        #print("\n\nVisiting Power: " + o.__repr__())

        # Get base and exponent.
        base, expo = o.operands()
        #print("\nbase: " + repr(base))
        #print("\nexponent: " + str(expo))

        # Visit base to get base code.
        base_code = self.visit(base)
        #print("base_code: " + str(base_code))

        # TODO: Are these safety checks needed? Need to check for None?
        if not () in base_code and len(base_code) != 1:
            error("Only support function type base: " + str(base_code))

        # Get the base code.
        val = base_code[()]

        # Handle different exponents
        if isinstance(expo, IntValue):
            return {(): self.format["power"](val, expo.value())}
        elif isinstance(expo, FloatValue):
            return {(): self.format["std power"](val, self.format["floating point"](expo.value()))}
        elif isinstance(expo, Function):
            exp = self.visit(expo)
            return {(): self.format["std power"](val, exp[()])}
        else:
            error("power does not support this exponent: " + repr(expo))

    def abs(self, o, *operands):
        #print("\n\nVisiting Abs: " + o.__repr__() + "with operands: " + "\n".join(map(str,operands)))

        # Prefetch formats to speed up code generation.
        format_abs = self.format["absolute value"]

        # TODO: Are these safety checks needed? Need to check for None?
        if len(operands) != 1 and not () in operands[0] and len(operands[0]) != 1:
            error("Abs expects one operand of function type: " + str(operands))

        # Take absolute value of operand.
        return {():format_abs(operands[0][()])}

    # -------------------------------------------------------------------------
    # Constant values (constantvalue.py).
    # -------------------------------------------------------------------------
    def format_scalar_value(self, value):
        #print("create_scalar_value: %d" % value)
        if value is None:
            return {():None}
        # TODO: Handle value < 0 better such that we don't have + -2 in the code.
        return {():self.format["floating point"](value)}

    # -------------------------------------------------------------------------
    # Constants (function.py).
    # -------------------------------------------------------------------------
    def create_constant_coefficient(self, count, component):
        coefficient = self.format["coeff"] + self.format["matrix access"](count, component)
        #print("create_constant_coefficient: " + coefficient)
        return {():coefficient}

    # -------------------------------------------------------------------------
    # FacetNormal (geometry.py).
    # -------------------------------------------------------------------------
    def facet_normal(self, o,  *operands):
        #print("Visiting FacetNormal:")

        # Get the component
        components = self.component()

        # Safety checks.
        if operands:
            error("Didn't expect any operands for FacetNormal: " + str(operands))
        if len(components) != 1:
            error("FacetNormal expects 1 component index: " + str(components))

        # We get one component.
        normal_component = self.format["normal component"](self.restriction, components[0])
        self.trans_set.add(normal_component)
        #print("Facet Normal Component: " + normal_component)
        return {():normal_component}

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
        #print("operand: " + str(operand))
        return operand

    # -------------------------------------------------------------------------
    # Helper functions for BasisFunction and Function).
    # -------------------------------------------------------------------------
    def __apply_transform(self, function, derivatives, multi):
        "Apply transformation (from derivatives) to basis or function."
        format_mult          = self.format["multiply"]
        format_transform     = self.format["transform"]

        # Add transformation if needed.
        transforms = []
        for i, direction in enumerate(derivatives):
            ref = multi[i]
            t = format_transform("JINV", ref, direction, self.restriction)
            self.trans_set.add(t)
            transforms.append(t)

        # Only multiply by basis if it is present.
        if function:
            prods = transforms + [function]
        else:
            prods = transforms

        return self.format["multiply"](prods)

    def create_basis_function(self, ufl_basis_function, derivatives, component, local_comp,
                  local_offset, ffc_element, transformation, multiindices):
        "Create code for basis functions, and update relevant tables of used basis."

        # Prefetch formats to speed up code generation.
        format_group         = self.format["grouping"]
        format_add           = self.format["add"]
        format_mult          = self.format["multiply"]
        format_transform     = self.format["transform"]
        format_detJ          = self.format["determinant"]
        format_inv           = self.format["inverse"]

        code = {}
        # Handle affine mappings.
        if transformation == AFFINE:
            # Loop derivatives and get multi indices.
            for multi in multiindices:
                deriv = [multi.count(i) for i in range(self.geo_dim)]
                if not any(deriv):
                    deriv = []
                # Call function to create mapping and basis name.
                mapping, basis = self.__create_mapping_basis(component, deriv, ufl_basis_function, ffc_element)
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
                    mapping, basis = self.__create_mapping_basis(c + local_offset, deriv, ufl_basis_function, ffc_element)
                    if basis is None:
                        if not mapping in code:
                            code[mapping] = []
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
                    if mapping in code:
                        code[mapping].append(self.__apply_transform(basis, derivatives, multi))
                    else:
                        code[mapping] = [self.__apply_transform(basis, derivatives, multi)]

        # Add sums and group if necessary.
        for key, val in code.items():
            if len(val) > 1:
                code[key] = format_group(format_add(val))
            elif val:
                code[key] = val[0]
            else:
                # Return a None (zero) because val == []
                code[key] = None

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

        basis = "1"
        if zeros and self.optimise_options["ignore zero tables"]:
            basis = None

        # If the loop index range is one we can look up the first component
        # in the psi array. If we only have ones we don't need the basis.
        if self.optimise_options["ignore ones"] and loop_index_range == 1 and ones:
            loop_index = "0"
        elif not basis is None:
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

    def create_function(self, ufl_function, derivatives, component, local_comp,
                  local_offset, ffc_element, quad_element, transformation, multiindices):
        "Create code for basis functions, and update relevant tables of used basis."

        # Prefetch formats to speed up code generation.
        format_mult          = self.format["multiply"]
        format_transform     = self.format["transform"]
        format_detJ          = self.format["determinant"]
        format_inv           = self.format["inverse"]

        code = []
        # Handle affine mappings.
        if transformation == AFFINE:
            # Loop derivatives and get multi indices.
            for multi in multiindices:
                deriv = [multi.count(i) for i in range(self.geo_dim)]
                if not any(deriv):
                    deriv = []
                # Call other function to create function name.
                function_name = self.__create_function_name(component, deriv, quad_element, ufl_function, ffc_element)
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
                    function_name = self.__create_function_name(c + local_offset, deriv, quad_element, ufl_function, ffc_element)
                    if function_name is None:
                        continue

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
                    code.append(self.__apply_transform(function_name, derivatives, multi))

        if not code:
            return None
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

        # If all basis are zero we just return None.
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

    # -------------------------------------------------------------------------
    # Helper functions for code_generation()
    # -------------------------------------------------------------------------
    def _count_operations(self, expression):
        return operation_count(expression, self.format)

    def _weight(self):
        # Create weight.
        weight = self.format["weight"](self.points)
        if self.points > 1:
            weight += self.format["array access"](self.format["integration points"])
        return weight

    def _create_value(self, val, weight, scale_factor):
        format_mult = self.format["multiply"]
        zero = False

        # Multiply value by weight and determinant
        value = format_mult([val, weight, scale_factor])

        return value, zero

