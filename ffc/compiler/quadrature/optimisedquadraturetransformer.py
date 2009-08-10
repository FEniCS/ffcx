"QuadratureTransformer (optimised) for quadrature code generation to translate UFL expressions."

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2009-03-18 -- 2009-08-08"
__copyright__ = "Copyright (C) 2009 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Python modules.
from numpy import shape, transpose

# UFL common.
from ufl.common import StackDict, product

# UFL Classes.
from ufl.classes import MultiIndex, FixedIndex

# UFL Algorithms.
from ufl.algorithms.transformations import Transformer
from ufl.algorithms import purge_list_tensors, expand_indices, propagate_restrictions
from ufl.algorithms.printing import tree_format

# FFC common modules.
from ffc.common.log import info, debug, error

# FFC compiler modules.
from ffc.compiler.tensor.multiindex import MultiIndex as FFCMultiIndex

# FFC fem modules.
from ffc.fem.createelement import create_element
from ffc.fem.finiteelement import AFFINE, CONTRAVARIANT_PIOLA, COVARIANT_PIOLA

# Utility and optimisation functions for quadraturegenerator.
from quadraturegenerator_utils import generate_loop, generate_psi_name, create_permutations
from quadraturetransformer import QuadratureTransformer
from symbolics import *

import time

class QuadratureTransformerOpt(QuadratureTransformer):
    "Transform UFL representation to quadrature code."

    def __init__(self, form_representation, domain_type, optimise_options, format):

        # Initialise base class.
        QuadratureTransformer.__init__(self, form_representation, domain_type, optimise_options, format)
        set_format(format)

    # -------------------------------------------------------------------------
    # Start handling UFL classes.
    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------
    # AlgebraOperators (algebra.py).
    # -------------------------------------------------------------------------
    def sum(self, o, *operands):
        debug("Visiting Sum: " + o.__repr__() + "\noperands: " + "\n".join(map(str, operands)))

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
            else:
                code[key] = val[0]
        return code

    def product(self, o, *operands):
        debug("\n\nVisiting Product: " + o.__repr__() + "with operands: " + "\n".join(map(str,operands)))

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
                if tuple(l) in code:
                    error("This key should not be in the code.")
                code[tuple(l)] = create_product(val + not_permute)
        else:
            code[()] = create_product(not_permute)

        return code

    def division(self, o, *operands):
        debug("\n\nVisiting Division: " + o.__repr__() + "with operands: " + "\n".join(map(str,operands)))

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
            debug("\nnum: " + repr(val))
            debug("\nden: " + repr(denominator))
            numerator_code[key] = create_fraction(val, denominator)
            debug("\ncode: " + str(numerator_code[key]))

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

        # Get the base code and create power.
        val = base_code.pop(())
        power = create_product([val]*expo.value())
        return {(): power}

    def abs(self, o, *operands):
        debug("\n\nVisiting Abs: " + o.__repr__() + "with operands: " + "\n".join(map(str,operands)))

        # TODO: Are these safety checks needed?
        if len(operands) != 1 and not () in operands[0] and len(operands[0]) != 1:
            error("Abs expects one operand of function type: " + str(operands))

        # Take absolute value of operand.
        operand = operands[0]
        for key, val in operand.items():
            new_val = create_symbol(self.format["absolute value"](str(val)), val.t)
            new_val.base_expr = val
            new_val.base_op = 1 # Add one operation for taking the absolute value.
            operand[key] = new_val
        return operand

    # -------------------------------------------------------------------------
    # Constant values (constantvalue.py).
    # -------------------------------------------------------------------------
    def scalar_value(self, o, *operands):
        "ScalarValue covers IntValue and FloatValue"
        debug("\n\nVisiting ScalarValue:" + o.__repr__())

        # FIXME: Might be needed because it can be IndexAnnotated?
        if operands:
            error("Did not expect any operands for ScalarValue: " + str((o, operands)))

        return {(): create_float(o.value())}

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
        return {(): create_symbol(coefficient, GEO)}

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
        return {(): create_symbol(coefficient, GEO)}

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
        return {(): create_symbol(coefficient, GEO)}

    # -------------------------------------------------------------------------
    # MathFunctions (mathfunctions.py).
    # -------------------------------------------------------------------------
    def __math_function(self, operands, format_function):
        # TODO: Are these safety checks needed?
        if len(operands) != 1 and not () in operands[0] and len(operands[0]) != 1:
            error("MathFunctions expect one operand of function type: " + str(operands))
        # Use format function on value of operand.
        operand = operands[0]
        for key, val in operand.items():
            new_val = create_symbol(format_function(str(val)), val.t)
            new_val.base_expr = val
            new_val.base_op = 1 # Add one operation for the math function.
            operand[key] = new_val
        debug("operand: " + str(operand))
        return operand

    def sqrt(self, o, *operands):
        debug("\n\nVisiting Sqrt: " + o.__repr__() + "with operands: " + "\n".join(map(str,operands)))
        # Call common math function.
        return self.__math_function(operands, self.format["sqrt"])

    def exp(self, o, *operands):
        debug("\n\nVisiting Exp: " + o.__repr__() + "with operands: " + "\n".join(map(str,operands)))
        # Call common math function.
        return self.__math_function(operands, self.format["exp"])

    def ln(self, o, *operands):
        debug("\n\nVisiting Ln: " + o.__repr__() + "with operands: " + "\n".join(map(str,operands)))
        # Call common math function.
        return self.__math_function(operands, self.format["ln"])

    def cos(self, o, *operands):
        debug("\n\nVisiting Cos: " + o.__repr__() + "with operands: " + "\n".join(map(str,operands)))
        # Call common math function.
        return self.__math_function(operands, self.format["cos"])

    def sin(self, o, *operands):
        debug("\n\nVisiting Sin: " + o.__repr__() + "with operands: " + "\n".join(map(str,operands)))
        # Call common math function.
        return self.__math_function(operands, self.format["sin"])

    # -------------------------------------------------------------------------
    # Helper functions for BasisFunction and Function).
    # -------------------------------------------------------------------------
    def create_basis_function(self, ufl_basis_function, component, derivatives):
        "Create code for basis functions, and update relevant tables of used basis."

        # Prefetch formats to speed up code generation.
        format_transform     = self.format["transform"]
        format_detJ          = self.format["determinant"]

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

            if ffc_element.component_element(component)[0].mapping() == AFFINE:
                # Call function to create mapping and basis name.
                mapping, basis = self.__create_mapping_basis(component, deriv, ufl_basis_function, ffc_element)

                # Add transformation if needed.
                transforms = []
                for i, direction in enumerate(derivatives):
                    ref = multi[i]
                    t = format_transform("JINV", ref, direction, self.restriction)
                    transforms.append(t)

                if mapping in code:
                    code[mapping].append(create_product([create_symbol(t, GEO) for t in transforms] + [basis]))
                else:
                    code[mapping] = [create_product([create_symbol(t, GEO) for t in transforms] + [basis])]
            # Handle non-affine mappings.
            else:
                for c in range(geo_dim):
                    # Call function to create mapping and basis name.
                    mapping, basis = self.__create_mapping_basis(c + local_offset, deriv, ufl_basis_function, ffc_element)

                    # Multiply basis by appropriate transform.
                    if ffc_element.component_element(component)[0].mapping() == COVARIANT_PIOLA:
                        dxdX = create_symbol(format_transform("JINV", c, local_comp, self.restriction), GEO)
                        basis = create_product([dxdX, basis])
                    elif ffc_element.component_element(component)[0].mapping() == CONTRAVARIANT_PIOLA:
                        detJ = create_fraction(create_float(1), create_symbol(format_detJ(self.restriction), GEO))
                        dXdx = create_symbol(format_transform("J", c, local_comp, self.restriction), GEO)
                        basis = create_product([detJ, dXdx, basis])
                    else:
                        error("Transformation is not supported: " + str(ffc_element.component_element(component)[0].mapping()))

                    # Add transformation if needed.
                    transforms = []
                    for i, direction in enumerate(derivatives):
                        ref = multi[i]
                        t = format_transform("JINV", ref, direction, self.restriction)
                        transforms.append(t)

                    if mapping in code:
                        code[mapping].append(create_product([create_symbol(t, GEO) for t in transforms] + [basis]))
                    else:
                        code[mapping] = [create_product([create_symbol(t, GEO) for t in transforms] + [basis])]

        # Add sums and group if necessary.
        for key, val in code.items():
            if len(val) > 1:
                code[key] = create_sum(val)
            else:
                code[key] = val[0]

        return code

    def __create_mapping_basis(self, component, deriv, ufl_basis_function, ffc_element):
        "Create basis name and mapping from given basis_info."

        # Get string for integration point.
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

        # Append the name to the set of used tables and create matrix access.
        basis = "0"
        if zeros and (self.optimise_options["ignore zero tables"] or self.optimise_options["remove zero terms"]):
            basis = create_float(0)
        elif self.optimise_options["ignore ones"] and loop_index_range == 1 and ones:
            basis = create_float(1)
            loop_index = "0"
        else:
            basis = create_symbol(name + basis_access, BASIS)
            self.psi_tables_map[basis] = name

        # Create the correct mapping of the basis function into the local element tensor.
        basis_map = loop_index
        if non_zeros and basis_map == "0":
            basis_map = str(non_zeros[1][0])
        elif non_zeros:
            basis_map = self.format["nonzero columns"](non_zeros[0]) + self.format["array access"](basis_map)
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
        format_transform     = self.format["transform"]
        format_detJ          = self.format["determinant"]

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
            if ffc_element.component_element(component)[0].mapping() == AFFINE:
                # Call other function to create function name.
                function_name = self.__create_function_name(component, deriv, quad_element, ufl_function, ffc_element)
                if not function_name:
                    continue

                # Add transformation if needed.
                transforms = []
                for i, direction in enumerate(derivatives):
                    ref = multi[i]
                    t = format_transform("JINV", ref, direction, self.restriction)
                    transforms.append(t)

                # Multiply function value by the transformations and add to code.
                code.append(create_product([create_symbol(t, GEO) for t in transforms] + [function_name]))

            # Handle non-affine mappings.
            else:
                for c in range(geo_dim):
                    function_name = self.__create_function_name(c + local_offset, deriv, quad_element, ufl_function, ffc_element)

                    # Multiply basis by appropriate transform.
                    if ffc_element.component_element(component)[0].mapping() == COVARIANT_PIOLA:
                        dxdX = create_symbol(format_transform("JINV", c, local_comp, self.restriction), GEO)
                        function_name = create_product([dxdX, function_name])
                    elif ffc_element.component_element(component)[0].mapping() == CONTRAVARIANT_PIOLA:
                        detJ = create_fraction(create_float(1), create_symbol(format_detJ(self.restriction), GEO))
                        dXdx = create_symbol(format_transform("J", c, local_comp, self.restriction), GEO)
                        function_name = create_product([detJ, dXdx, function_name])
                    else:
                        error("Transformation is not supported: ", str(ffc_element.component_element(component)[0].mapping()))

                    # Add transformation if needed.
                    transforms = []
                    for i, direction in enumerate(derivatives):
                        ref = multi[i]
                        t = format_transform("JINV", ref, direction, self.restriction)
                        self.trans_set.add(t)
                        transforms.append(t)

                    # Multiply function value by the transformations and add to code.
                    code.append(create_product([create_symbol(t, GEO) for t in transforms] + [function_name]))
        if not code:
            return "0"
        elif len(code) > 1:
            code = create_sum(code)
        else:
            code = code[0]

        return code

    def __create_function_name(self, component, deriv, quad_element, ufl_function, ffc_element):

        # Prefetch formats to speed up code generation.
        format_ip            = self.format["integration points"]

        # Pick first free index of secondary type
        # (could use primary indices, but it's better to avoid confusion).
        loop_index = self.format["free secondary indices"][0]

        # Create basis access, we never need to map the entry in the basis table
        # since we will either loop the entire space dimension or the non-zeros.
        if self.points == 1:
            format_ip = "0"
        basis_access = self.format["matrix access"](format_ip, loop_index)
        ACCESS = IP

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
        if zeros and self.optimise_options["ignore zero tables"]:
            return create_float(0)

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
            # TODO: place this in basis map first and the use only add to map
            # if the function value is actually used.
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
            ACCESS = GEO
        except:
            pass

        coefficient = create_symbol(self.format["coeff"] +\
                      self.format["matrix access"](str(ufl_function.count()), coefficient_access), ACCESS)
        function_expr = coefficient
        if basis_name:
            function_expr = create_product([create_symbol(basis_name, ACCESS), coefficient])

        # If we have a quadrature element (or if basis was deleted) we don't need the basis.
        if quad_element or not basis_name:
            function_name = coefficient
        else:
            # Check if the expression to compute the function value is already in
            # the dictionary of used function. If not, generate a new name and add.
            function_name = create_symbol(self.format["function value"] + str(self.function_count), ACCESS)
            if not function_expr in self.functions:
                self.functions[function_expr] = (function_name, loop_index_range)
                # Increase count.
                self.function_count += 1
            else:
                function_name, index_r = self.functions[function_expr]
                # Check just to make sure.
                if not index_r == loop_index_range:
                    error("Index ranges does not match")
        return function_name

def generate_code(integrand, transformer, Indent, format, interior):
    """Generate code from a UFL integral type. It generates all the code that
    goes inside the quadrature loop."""

    # Prefetch formats to speed up code generation.
    format_comment          = format["comment"]
    format_float_decl       = format["float declaration"]
    format_F                = format["function value"]
    format_float            = format["floating point"]
    format_add_equal        = format["add equal"]
    format_nzc              = format["nonzero columns"](0).split("0")[0]
    format_r                = format["free secondary indices"][0]
    format_array_access     = format["array access"]
    format_scale_factor     = format["scale factor"]
    format_add              = format["add"]
    format_mult             = format["multiply"]
    format_tensor           = format["element tensor quad"]
    format_Gip              = format["geometry tensor"] + format["integration points"]

    # Initialise return values.
    code = []
    num_ops = 0

    debug("\nQG, Using Transformer.")

    # Apply basic expansions.
    # TODO: Figure out if there is a 'correct' order of doing this.
    # In form.form_data().form, which we should be using, coefficients have
    # been mapped and derivatives expanded. So it should be enough to just
    # expand_indices and purge_list_tensors.
    new_integrand = expand_indices(integrand)
    new_integrand = purge_list_tensors(new_integrand)
    # Only propagate restrictions if we have an interior integral.
    if interior:
        new_integrand = propagate_restrictions(new_integrand)
    debug("\nExpanded integrand\n" + str(tree_format(new_integrand)))
    # Let the Transformer create the loop code.
    info("Transforming UFL integrand...")
    t = time.time()
    loop_code = transformer.visit(new_integrand)
    info("done, time = %f" % (time.time() - t))

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
            number = int(str(name).strip(format_F))
            # TODO: This check can be removed for speed later.
            if number in function_numbers:
                error("This is definitely not supposed to happen!")
            function_numbers.append(number)
            # Get number of operations to compute entry and add to function operations count.
            f_ops = function.ops() + 1
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
    ACCESS = GEO
    weight = format["weight"](transformer.points)
    if transformer.points > 1:
        weight += format_array_access(format["integration points"])
        ACCESS = IP

    # Generate entries, multiply by weights and sort after primary loops.
    loops = {}
    ip_consts = {}

    info("Optimising code...")
    t = time.time()
    for key, val in loop_code.items():
        # If value was zero continue.
        if val == None:
            continue
        # Multiply by weight and determinant, add both to set of used weights and transforms.
        value = create_product([val, create_symbol(weight, ACCESS), create_symbol(format_scale_factor, GEO)])
        value = optimise_code(value, ip_consts, transformer.geo_consts, transformer.trans_set)

        # Only continue if value is not zero.
        if not value.val:
            continue

        # Add weight and determinant to sets.
        transformer.used_weights.add(transformer.points)
        transformer.trans_set.add(format_scale_factor)

        # Update the set of used psi tables through the name map.
        transformer.used_psi_tables.update([transformer.psi_tables_map[b] for b in value.get_unique_vars(BASIS)])

        # Compute number of operations to compute entry and create comment
        # (add 1 because of += in assignment).
        entry_ops = value.ops() + 1
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
    info("           done, time = %f" % (time.time() - t))

    info("Writing code...")
    t = time.time()
    # Generate code for ip constant declarations.
    ip_const_ops, ip_const_code = generate_aux_constants(ip_consts, format_Gip,\
                                    format["const float declaration"], True)
    num_ops += ip_const_ops
    if ip_const_code:
        code += ["", format["comment"]("Number of operations to compute ip constants: %d" %ip_const_ops)]
        code += ip_const_code

    # Write all the loops of basis functions.
    for loop, ops_lines in loops.items():
        ops, lines = ops_lines

        # Add number of operations for current loop to total count.
        num_ops += ops
        code += ["", format_comment("Number of operations for primary indices = %d" % ops)]
        code += generate_loop(lines, loop, Indent, format)
    info("              done, time = %f" % (time.time() - t))

    return (code, num_ops)


