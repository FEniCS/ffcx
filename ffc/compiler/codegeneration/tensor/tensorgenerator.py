"Code generator for tensor representation"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-03 -- 2008-06-12"
__copyright__ = "Copyright (C) 2004-2008 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian B. Oelgaard 2007
# Modified by Marie Rognes (meg@math.uio.no) 2007

# Python modules
from sets import Set

# FFC common modules
from ffc.common.constants import *

# FFC language modules
from ffc.compiler.language.index import *

# FFC code generation common modules
from ffc.compiler.codegeneration.common.codegenerator import *

# FFC format modules
from ffc.compiler.format.removeunused import *

class TensorGenerator(CodeGenerator):
    "Code generator for for tensor representation"

    def __init__(self):
        "Constructor"

        # Initialize common code generator
        CodeGenerator.__init__(self)

    def generate_cell_integral(self, form_representation, sub_domain, format):
        """Generate dictionary of code for cell integral from the given
        form representation according to the given format"""

        # Extract terms for sub domain
        terms = [term for term in form_representation.cell_tensor if term.monomial.integral.sub_domain == sub_domain]

        # Special case: zero contribution
        if len(terms) == 0:
            element_code = self.__generate_zero_element_tensor(form_representation.cell_tensor, format)
            return {"tabulate_tensor": element_code, "members": ""}

        debug("")
        # Generate element code + set of used geometry terms
        element_code, geo_set, tensor_ops = self.__generate_element_tensor(terms, format)

        # Generate geometry code + set of used coefficients + set of jacobi terms
        geo_code, coeff_set, trans_set, geo_ops = self.__generate_geometry_tensors(terms, geo_set, format)
        total_ops = tensor_ops + geo_ops

        # Generate code for manipulating coefficients
        coeff_code = self.__generate_coefficients(terms, coeff_set, format) 

        # Get Jacobian snippet
        jacobi_code = [format["generate jacobian"](form_representation.cell_dimension, Integral.CELL)]

        # Remove unused declarations
        code = self.__remove_unused(jacobi_code, trans_set, format)

        # Add coefficient and geometry tensor declarations
        code.append(format["comment"]("Number of operations to compute element tensor = %d" % total_ops))
        code += coeff_code + geo_code

        # Add element code
        code += [""] + [format["comment"]("Compute element tensor")]
        code += element_code

        debug("Number of operations to compute tensor: %d" % total_ops)

        return {"tabulate_tensor": code, "members": ""}

    def generate_exterior_facet_integral(self, form_representation, sub_domain, format):
        """Generate dictionary of code for exterior facet integral from the given
        form representation according to the given format"""

        # Extract terms for sub domain
        terms = [[term for term in t if term.monomial.integral.sub_domain == sub_domain] for t in form_representation.exterior_facet_tensors]

        # Special case: zero contribution
        if all([len(t) == 0 for t in terms]):
            element_code = self.__generate_zero_element_tensor(form_representation.exterior_facet_tensors[0], format)
            return {"tabulate_tensor": (element_code, []), "members": ""}

        num_facets = len(terms)
        cases = [None for i in range(num_facets)]

        # Generate element code + set of used geometry terms
        geo_set = Set()
        debug("")
        tensor_ops = 0
        for i in range(num_facets):
            case, g_set, tensor_ops = self.__generate_element_tensor(terms[i], format)
            cases[i] = case
            geo_set = geo_set | g_set
            debug("Number of operations to compute element tensor for facet %d: %d"% (i, tensor_ops))

        # Generate code for geometry tensor (should be the same so pick first)
        # Generate set of used coefficients + set of jacobi terms
        geo_code, coeff_set, trans_set, geo_ops = self.__generate_geometry_tensors(terms[0], geo_set, format)
        debug("Number of operations to compute geometry terms (should be added): %d" % geo_ops)
        total_ops = tensor_ops + geo_ops

        # Generate code for manipulating coefficients (should be the same so pick first)
        coeff_code = self.__generate_coefficients(terms[0], coeff_set, format)

        # Get Jacobian snippet
        jacobi_code = [format["generate jacobian"](form_representation.cell_dimension, Integral.EXTERIOR_FACET)]

        # Remove unused declarations
        code = self.__remove_unused(jacobi_code, trans_set, format)

        code.append(format["comment"]("Number of operations to compute element tensor = %d" % total_ops))

        # Add coefficient and geometry tensor declarations
        code += coeff_code + geo_code

        # Add element code
        code += [""] + [format["comment"]("Compute element tensor for all facets")]

        return {"tabulate_tensor": (code, cases), "members": ""}
    
    def generate_interior_facet_integral(self, form_representation, sub_domain, format):
        """Generate dictionary of code for interior facet integral from the given
        form representation according to the given format"""

        # Extract terms for sub domain
        terms = [[[term for term in t2 if term.monomial.integral.sub_domain == sub_domain] for t2 in t1] for t1 in form_representation.interior_facet_tensors]

        # Special case: zero contribution
        if all([len(t) == 0 for tt in terms for t in tt]):
            element_code = self.__generate_zero_element_tensor(form_representation.interior_facet_tensors[0][0], format)
            return {"tabulate_tensor": (element_code, []), "members": ""}

        num_facets = len(terms)
        cases = [[None for j in range(num_facets)] for i in range(num_facets)]

        # Generate element code + set of used geometry terms
        geo_set = Set()
        debug("")
        tensor_ops = 0
        for i in range(num_facets):
            for j in range(num_facets):
                case, g_set, tensor_ops = self.__generate_element_tensor(terms[i][j], format)
                cases[i][j] = case
                geo_set = geo_set | g_set
                debug("Number of operations to compute element tensor for facets (%d, %d): %d" % (i, j, tensor_ops))

        # Generate code for geometry tensor (should be the same so pick first)
        # Generate set of used coefficients + set of jacobi terms
        geo_code, coeff_set, trans_set, geo_ops = self.__generate_geometry_tensors(terms[0][0], geo_set, format)
        debug("Number of operations to compute geometry terms (should be added): %d" % geo_ops)
        total_ops = tensor_ops + geo_ops

        # Generate code for manipulating coefficients (should be the same so pick first)
        coeff_code = self.__generate_coefficients(terms[0][0], coeff_set, format)

        # Get Jacobian snippet
        jacobi_code = [format["generate jacobian"](form_representation.cell_dimension, Integral.INTERIOR_FACET)]

        # Remove unused declarations
        code = self.__remove_unused(jacobi_code, trans_set, format)

        code.append(format["comment"]("Number of operations to compute element tensor = %d" % total_ops))

        # Add coefficient and geometry tensor declarations
        code += coeff_code + geo_code

        # Add element code
        code += [""] + [format["comment"]("Compute element tensor for all facet-facet combinations")]

        return {"tabulate_tensor": (code, cases), "members": ""}

    def __generate_coefficients(self, terms, coeff_set, format):
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
                    index = coefficient.n0.index
                    if term.monomial.integral.type == Integral.INTERIOR_FACET:
                        space_dimension = 2*len(coefficient.index.range)
                    else:
                        space_dimension = len(coefficient.index.range)
                    for l in range(space_dimension):
                        # If coefficient is not used don't declare it
                        if not format["modified coefficient access"](index, j, k, l) in coeff_set:
                            continue
                        name = format["modified coefficient declaration"](index, j, k, l)
                        value = format["coefficient"](index, l)
                        for l in range(len(coefficient.ops)):
                            op = coefficient.ops[len(coefficient.ops) - 1 - l]
                            if op == Operators.INVERSE:
                                value = format["inverse"](value)
                            elif op == Operators.MODULUS:
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

    def __generate_geometry_tensors(self, terms, geo_set, format):
        "Generate list of declarations for computation of geometry tensors"

        # Generate code as a list of declarations
        code = []    

        # Iterate over all terms
        j = 0
        coeff_set = Set()
        trans_set = Set()
        num_ops = 0
        for i in range(len(terms)):

            term = terms[i]

            # Get list of secondary indices (should be the same so pick first)
            aindices = terms[i].G[0].a.indices

            # Iterate over secondary indices
            for a in aindices:

                # Skip code generation if term is not used
                if not format["geometry tensor access"](i,a) in geo_set:
                    continue

                # Compute factorized values
                values = []
                jj = j
                for G in term.G:
                    val, c_set, t_set, entry_ops = self.__generate_entry(G, a, jj, format)                    
                    values += [val]
                    num_ops += entry_ops
                    coeff_set = coeff_set | c_set
                    trans_set = trans_set | t_set
                    jj += 1

                # Sum factorized values
                if values:
                    num_ops += len(values) - 1
                name = format["geometry tensor declaration"](i, a)
                value = format["add"](values)

                # Multiply with determinant factor
                # FIXME: dets = pick_first([G.determinants for G in term.G])                
                dets = term.G[0].determinants
                value = self.__multiply_value_by_det(value, dets, format, len(values) > 1)
                num_ops += 1

                # Add determinant to transformation set
                if dets:
                    d0 = [format["power"](format["determinant"](det.restriction),
                                          det.power) for det in dets]
                    trans_set.add(format["multiply"](d0))

                # Add declaration
                code += [(name, value)]

            j += len(term.G)

        # Add comments
        code = [format["comment"]("Compute geometry tensors"), format["comment"]("Number of operations to compute decalrations = %d" %num_ops)] + code

        # Add scale factor
        trans_set.add(format["scale factor"])

        return (code, coeff_set, trans_set, num_ops)

    def __generate_element_tensor(self, terms, format):
        "Generate list of declarations for computation of element tensor"

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
        geo_set = Set()
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
                            geo_set.add(gk)
                            num_ops += 1
                        elif value:
                            value = format_add([value, format_multiply([format_floating_point(a0), gk])])
                            geo_set.add(gk)
                            num_ops += 1
                        else:
                            value = format_multiply([format_floating_point(a0), gk])
                            geo_set.add(gk)
                        num_ops += 1
                    else:
                        num_dropped += 1
            value = value or zero
            code += [(name, value)]
            k += 1

        print geo_set

        code = [format["comment"]("Number of operations to compute tensor = %d" %num_ops)] + code
        return (code, geo_set, num_ops)

    def __generate_zero_element_tensor(self, A, format):
        "Generate list of declarations for zero element tensor"

        # Generate code as a list of declarations
        code = []    
    
        # Get indices (check first term)
        indices = A[0].A0.i.indices

        # Prefetch formats to speed up code generation
        format_element_tensor  = format["element tensor"]
        format_floating_point  = format["floating point"]

        # Set entries to zero
        for k in range(len(indices)):
            name = format_element_tensor(indices[k], k)
            value = format_floating_point(0.0)
            code += [(name, value)]

        return code

    def __generate_entry(self, G, a, i, format):
        "Generate code for the value of entry a of geometry tensor G"

        coeff_set = Set()
        trans_set = Set()

        # Compute product of factors outside sum
        factors = []
        num_ops = 0
        for j in range(len(G.coefficients)):
            c = G.coefficients[j]
            if not c.index.type == Index.AUXILIARY_G:
                coefficient = format["modified coefficient access"](c.n1.index, i, j, c.index([], a, [], []))
                coeff_set.add(coefficient)
                factors += [coefficient]
            
        for t in G.transforms:
            if not (t.index0.type == Index.AUXILIARY_G or t.index1.type == Index.AUXILIARY_G):
                trans = format["transform"](t.type, t.index0([], a, [], []), \
                                            t.index1([], a, [], []), \
                                            t.restriction)
                factors += [trans]
                trans_set.add(trans)
        if factors:
            num_ops += len(factors) - 1

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
                    coeff_set.add(coefficient)
                    factors += [coefficient]
            for t in G.transforms:
                if t.index0.type == Index.AUXILIARY_G or t.index1.type == Index.AUXILIARY_G:
                    trans = format["transform"](t.type, t.index0([], a, [], b), \
                                                            t.index1([], a, [], b), \
                                                            t.restriction)
                    factors += [trans]
                    trans_set.add(trans)
            if factors:
                num_ops += len(factors) - 1
            terms += [format["multiply"](factors)]

        if terms:
            num_ops += len(terms) - 1

        sum = format["add"](terms)
        if sum: sum = format["grouping"](sum)
        if sum: f1 = [sum]
        else: f1 = []

        fs = f0 + f1
        if not fs: fs = ["1.0"]
        else:
            num_ops += len(fs) - 1

        # Compute product of all factors
        return (format["multiply"](fs), coeff_set, trans_set, num_ops)

    def __multiply_value_by_det(self, value, dets, format, is_sum):
        if dets:
            d0 = [format["power"](format["determinant"](det.restriction),
                                  det.power) for det in dets]
        else:
            d0 = []
        if value == "1.0":
            v = []
        elif is_sum:
            v = [format["grouping"](value)]
        else:
            v = [value]
        return format["multiply"](d0 + [format["scale factor"]] + v)

    def __remove_unused(self, code, set, format):
        "Remove unused variables so that the compiler will not complain"

        # Normally, the removal of unused variables should happen at the
        # formatting stage, but since the code for the tensor contraction
        # may grow to considerable size, we make an exception and remove
        # unused variables here when we know the names of the used
        # variables. No searching necessary and much, much, much faster.
        
        if code:
            # Generate body of code, using the format
            lines = format["generate body"](code)

            # Generate auxiliary code line that uses all members of the set (to trick remove_unused)
            line_set = format["add equal"]("A", format["multiply"](set))
            lines += "\n" + line_set

            # Remove unused Jacobi declarations
            code = remove_unused(lines)

            # Delete auxiliary line
            code = code.replace("\n" + line_set, "")

            return [code]
        else:
            return code
