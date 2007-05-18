"Code generator for tensor representation"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-03 -- 2007-05-10"
__copyright__ = "Copyright (C) 2004-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC common modules
from ffc.common.constants import *

# FFC language modules
from ffc.compiler.language.index import *

# FFC code generation common modules
from ffc.compiler.codegeneration.common.codegenerator import *

class TensorGenerator(CodeGenerator):
    "Code generator for for tensor representation"

    def __init__(self):
        "Constructor"

        # Initialize common code generator
        CodeGenerator.__init__(self)

    def generate_cell_integral(self, form_representation, sub_domain, format):
        """Generate dictionary of code for cell integral from the given
        form representation according to the given format"""

        # Extract terms
        terms = form_representation.cell_tensor
        if len(terms) == 0:
            return None

        # Generate code for manipulating coefficients
        code = self.__generate_coefficients(terms, format)

        # Generate code for geometry tensor
        code += self.__generate_geometry_tensors(terms, format)

        # Generate code for sign changes
        (sign_code, change_signs) = self.__generate_signs(terms, format)
        code += sign_code
        
        # Generate code for element tensor(s)
        code += [""] + [format["comment"]("Compute element tensor")]
        code += self.__generate_element_tensor(terms, change_signs, format)

        return {"tabulate_tensor": code}

    def generate_exterior_facet_integral(self, form_representation, sub_domain, format):
        """Generate dictionary of code for exterior facet integral from the given
        form representation according to the given format"""

        # Extract terms
        terms = form_representation.exterior_facet_tensors
        if len(terms) == 0:
            return None

        # Generate code for manipulating coefficients (should be the same so pick first)
        code = self.__generate_coefficients(terms[0], format)
        
        # Generate code for geometry tensor (should be the same so pick first)
        code += self.__generate_geometry_tensors(terms[0], format)

        # Generate code for element tensor(s)
        code += [""] + [format["comment"]("Compute element tensor for all facets")]
        num_facets = len(terms)
        cases = [None for i in range(num_facets)]
        for i in range(num_facets):
            cases[i] = self.__generate_element_tensor(terms[i], False, format)

        return {"tabulate_tensor": (code, cases)}
    
    def generate_interior_facet_integral(self, form_representation, sub_domain, format):
        """Generate dictionary of code for interior facet integral from the given
        form representation according to the given format"""

        # Extract terms
        terms = form_representation.interior_facet_tensors
        if len(terms) == 0:
            return None

        # Generate code for manipulating coefficients (should be the same so pick first)
        code = self.__generate_coefficients(terms[0][0], format)
        
        # Generate code for geometry tensor (should be the same so pick first)
        code += self.__generate_geometry_tensors(terms[0][0], format)

        # Generate code for element tensor(s)
        code += [""] + [format["comment"]("Compute element tensor for all facet-facet combinations")]
        num_facets = len(terms)
        cases = [[None for j in range(num_facets)] for i in range(num_facets)]
        for i in range(num_facets):
            for j in range(num_facets):
                cases[i][j] = self.__generate_element_tensor(terms[i][j], False, format)

        return {"tabulate_tensor": (code, cases)}

    def __generate_coefficients(self, terms, format):
        "Generate code for manipulating coefficients"

        # Generate code as a list of declarations
        code = []

        # Add comment
        code += [format["comment"]("Compute coefficients")]

        # A coefficient is identified by 4 numbers:
        #
        #   0 - the number of the function
        #   1 - the position of the (factored) monomial it appears in
        #   2 - the position of the coefficient inside the monomial
        #   3 - the position of the expansion coefficient

        # Iterate over all terms
        j = 0

        for term in terms:
            for G in term.G:
                for k in range(len(G.coefficients)):
                    coefficient = G.coefficients[k]
                    if term.monomial.integral.type == Integral.INTERIOR_FACET:
                        space_dimension = 2*len(coefficient.n0.range)
                    else:
                        space_dimension = len(coefficient.n0.range)
                    for l in range(space_dimension):
                        name = format["modified coefficient declaration"](coefficient.n0.index, j, k, l)
                        value = format["coefficient"](coefficient.n0.index, l)
                        for l in range(len(coefficient.ops)):
                            op = coefficient.ops[len(coefficient.ops) - 1 - l]
                            if op == Operators.INVERSE:
                                value = format["inverse"](value)
                            elif op == Operators.ABS:
                                value = format["absolute value"](value)
                            elif op == Operators.SQRT:
                                value = format["sqrt"](value)
                        code += [(name, value)]
                j += 1

        # Don't add code if there are no coefficients
        if len(code) == 1:
            return []

        # Add newline
        code += [""]

        return code

    def __generate_geometry_tensors(self, terms, format):
        "Generate list of declarations for computation of geometry tensors"

        # Generate code as a list of declarations
        code = []    
        
        # Add comment
        code += [format["comment"]("Compute geometry tensors")]

        # Iterate over all terms
        j = 0
        for i in range(len(terms)):

            term = terms[i]

            # Get list of secondary indices (should be the same so pick first)
            aindices = terms[i].G[0].a.indices

            # Iterate over secondary indices
            for a in aindices:

                # Compute factorized values
                values = []
                jj = j
                for G in term.G:
                    values += [self.__generate_entry(G, a, jj, format)]
                    jj += 1

                # Sum factorized values
                name = format["geometry tensor declaration"](i, a)
                value = format["add"](values)

                # Multiply with determinant factor
                det = pick_first([G.determinant for G in term.G])
                value = self.__multiply_value_by_det(value, det, format, len(values) > 1)
                
                # Add declaration
                code += [(name, value)]

            j += len(term.G)

        return code
    
    def __generate_element_tensor(self, terms, sign_changes, format):
        "Generate list of declaration for computation of element tensor"

        # Generate code as a list of declarations
        code = []    
    
        # Get list of primary indices (should be the same so pick first)
        iindices = terms[0].A0.i.indices
    
        # Prefetch formats to speed up code generation
        format_element_tensor  = format["element tensor"]
        format_geometry_tensor = format["geometry tensor access"]
        format_add             = format["add"]
        format_subtract        = format["subtract"]
        format_multiply        = format["multiply"]
        format_floating_point  = format["floating point"]
        format_epsilon         = format["epsilon"]

        # Generate code for geometry tensor entries
        gk_tensor = [ ( [(format_geometry_tensor(j, a), a) for a in terms[j].A0.a.indices], j) for j in range(len(terms)) ]

        # Generate code for computing the element tensor
        k = 0
        num_dropped = 0
        num_ops = 0
        zero = format_floating_point(0.0)
        for i in iindices:
            name = format_element_tensor(i, k)
            value = None
            for (gka, j) in gk_tensor:
                A0 = terms[j].A0
                for (gk, a) in gka:
                    a0 = A0.A0[tuple(i + a)]
                    if abs(a0) > format_epsilon:
                        if value and a0 < 0.0:
                            value = format_subtract([value, format_multiply([format_floating_point(-a0), gk])])
                        elif value:
                            value = format_add([value, format_multiply([format_floating_point(a0), gk])])
                        else:
                            value = format_multiply([format_floating_point(a0), gk])
                        num_ops += 1
                    else:
                        num_dropped += 1

            # Add sign changes as appropriate.
            if sign_changes:
                value = self.__add_sign(value, 0, i, format)
            value = value or zero
            code += [(name, value)]
            k += 1

        return code
    
    def __generate_entry(self, G, a, i, format):
        "Generate code for the value of entry a of geometry tensor G"
    
        # Compute product of factors outside sum
        factors = []
        for j in range(len(G.coefficients)):
            c = G.coefficients[j]
            if not c.index.type == Index.AUXILIARY_G:
                coefficient = format["modified coefficient access"](c.n1.index, i, j, c.index([], a, [], []))
                factors += [coefficient]
        for t in G.transforms:
            if not (t.index0.type == Index.AUXILIARY_G or  t.index1.type == Index.AUXILIARY_G):
                factors += [format["transform"](t.type, t.index0([], a, [], []), \
                                                        t.index1([], a, [], []), \
                                                        t.restriction),]
        monomial = format["multiply"](factors)
        if monomial: f0 = [monomial]
        else: f0 = []
    
        # Compute sum of monomials inside sum
        terms = []
        for b in G.b.indices:
            factors = []
            for j in range(len(G.coefficients)):
                c = G.coefficients[j]
                if c.index.type == Index.AUXILIARY_G:
                    coefficient = format["modified coefficient access"](c.n1.index, i, j, c.index([], a, [], b))
                    factors += [coefficient]
            for t in G.transforms:
                if t.index0.type == Index.AUXILIARY_G or t.index1.type == Index.AUXILIARY_G:
                    factors += [format["transform"](t.type, t.index0([], a, [], b), \
                                                            t.index1([], a, [], b), \
                                                            t.restriction)]
            terms += [format["multiply"](factors)]
        sum = format["add"](terms)
        if sum: sum = format["grouping"](sum)
        if sum: f1 = [sum]
        else: f1 = []

        fs = f0 + f1
        if not fs: fs = ["1.0"]

        # Compute product of all factors
        return format["multiply"](fs)

    def __generate_signs(self, terms, format):
        "Generate list of declarations for computation of signs"
        code = []
        computed = {}
        for j in range(len(terms)):
            monomial = terms[j].monomial
            # Inspect each basis function (identified by its index)
            # and check whether sign changes are relevant.
            for basisfunction in monomial.basisfunctions:
                index = basisfunction.index
                if not str(index) in computed:
                    necessary = False
                    element = basisfunction.element
                    declarations = []
                    dof_entities = DofMap(element).dof_entities();

                    # Go through the topological entities associated
                    # with each basis function/dof. If the element is
                    # a piola mapped element and the basis function is
                    # associated with an edge, we calculate the
                    # possible sign change.
                    for no in dof_entities:
                        (entity, entity_no) = dof_entities[no]
                        name = format["sign tensor"](j, index.index, no)
                        if entity == 1 and element.space_mapping(no) == Mapping.PIOLA:
                            necessary = True
                            value = format["call edge sign"](entity_no)
                            # If the sign of this edge already has
                            # been computed, refer to that entry instead.
                            if value in computed:
                                value = computed[value]
                            else:
                                computed[value] = name
                        else:
                            value = "1"

                        # Add to declarations
                        declarations += [(format["sign tensor declaration"](name), value)]    

                    # Add declarations for this basis function to the code
                    code += declarations
                    computed[str(index)] = True
                    
        if necessary:
            code.insert(0, format["comment"]("Compute signs"))
            code.insert(0, format["snippet edge signs"](2))
            return (code, True)
        else:
            return ([], False) # Return [] is the case of no sign changes...)

    def __add_sign(self, value, j, i, format):
        if value:
            value = format["grouping"](value)
            for k in range(len(i)):
                value = format["multiply"]([format["sign tensor"](j, k, i[k]), value])
        return value

    def __multiply_value_by_det(self, value, det, format, is_sum):
        if det: d0 = [format["power"](format["determinant"], det)]
        else: d0 = []
        if value == "1.0":
            v = []
        elif is_sum:
            v = [format["grouping"](value)]
        else:
            v = [value]
        return format["multiply"](d0 + [format["scale factor"]] + v)
