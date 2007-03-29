"Code generator for tensor representation"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-03 -- 2007-03-30"
__copyright__ = "Copyright (C) 2004-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC common modules
from ffc.common.constants import *

# FFC language modules
from ffc.compiler.language.index import *

# FFC code generation modules
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
        
        # Generate code for geometry tensor
        code = self.__generate_geometry_tensors(terms, format)
        
        # Generate code for element tensor(s)
        code += [""] + [format["comment"]("Compute element tensor")]
        code += self.__generate_element_tensor(terms, format)

        return {"tabulate_tensor": code}

    def generate_exterior_facet_integral(self, form_representation, sub_domain, format):
        """Generate dictionary of code for exterior facet integral from the given
        form representation according to the given format"""

        # Extract terms
        terms = form_representation.exterior_facet_tensors
        if len(terms) == 0:
            return None
        
        # Generate code for geometry tensor (should be the same so pick first)
        common = self.__generate_geometry_tensors(terms[0], format)

        # Generate code for element tensor(s)
        common += [""] + [format["comment"]("Compute element tensor for all facets")]
        num_facets = len(terms)
        cases = [None for i in range(num_facets)]
        for i in range(num_facets):
            cases[i] = self.__generate_element_tensor(terms[i], format)

        return {"tabulate_tensor": (common, cases)}
    
    def generate_interior_facet_integral(self, form_representation, sub_domain, format):
        """Generate dictionary of code for interior facet integral from the given
        form representation according to the given format"""

        # Extract terms
        terms = form_representation.interior_facet_tensors
        if len(terms) == 0:
            return None
        
        # Generate code for geometry tensor (should be the same so pick first)
        common = self.__generate_geometry_tensors(terms[0][0], format)

        # Generate code for element tensor(s)
        common += [""] + [format["comment"]("Compute element tensor for all facet-facet combinations")]
        num_facets = len(terms)
        cases = [[None for j in range(num_facets)] for i in range(num_facets)]
        for i in range(num_facets):
            for j in range(num_facets):
                cases[i][j] = self.__generate_element_tensor(terms[i][j], format)

        return {"tabulate_tensor": (common, cases)}

    def __generate_geometry_tensors(self, terms, format):
        "Generate list of declarations for computation of geometry tensors"

        # Generate code as a list of declarations
        code = []    
        
        # Add comment
        code += [format["comment"]("Compute geometry tensors")]

        # Iterate over all terms
        for j in range(len(terms)):

            # Get list of secondary indices (should be the same so pick first)
            aindices = terms[j].G[0].a.indices

            # Iterate over secondary indices
            for a in aindices:

                # Sum factorized values
                name = format["geometry tensor declaration"](j, a)
                value = format["add"]([self.__generate_entry(G, a, format) for G in terms[j].G])

                # Add declaration
                code += [(name, value)]

        return code
    
    def __generate_element_tensor(self, terms, format):
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
                    if abs(a0) > FFC_EPSILON:
                        if value and a0 < 0.0:
                            value = format_subtract([value, format_multiply([format_floating_point(-a0), gk])])
                        elif value:
                            value = format_add([value, format_multiply([format_floating_point(a0), gk])])
                        else:
                            value = format_multiply([format_floating_point(a0), gk])
                        num_ops += 1
                    else:
                        num_dropped += 1
            value = value or zero
            code += [(name, value)]
            k += 1

        return code
    
    def __generate_entry(self, G, a, format):
        "Generate code for the value of entry a of geometry tensor G"
    
        # Compute product of factors outside sum
        factors = []
        for c in G.constants:
            if c.inverted:
                factors += ["(1.0/" + format["constant"](c.number.index) + ")"]
            else:
                factors += [format["constant"](c.number.index)]
        for c in G.coefficients:
            if not c.index.type == Index.AUXILIARY_G:
                coefficient = format["coefficient"](c.n1.index, c.index([], a, [], []))
                factors += [coefficient]
        for t in G.transforms:
            if not (t.index0.type == Index.AUXILIARY_G or  t.index1.type == Index.AUXILIARY_G):
                factors += [format["inverse transform"](t.index0([], a, [], []), \
                                                        t.index1([], a, [], []), \
                                                        t.restriction),]
        monomial = format["multiply"](factors)
        if monomial: f0 = [monomial]
        else: f0 = []
    
        # Compute sum of monomials inside sum
        terms = []
        for b in G.b.indices:
            factors = []
            for c in G.coefficients:
                if c.index.type == Index.AUXILIARY_G:
                    coefficient = format["coefficient"](c.n1.index, c.index([], a, [], b))
                    factors += [coefficient]
            for t in G.transforms:
                if t.index0.type == Index.AUXILIARY_G or t.index1.type == Index.AUXILIARY_G:
                    factors += [format["inverse transform"](t.index0([], a, [], b), \
                                                            t.index1([], a, [], b), \
                                                            t.restriction)]
            terms += [format["multiply"](factors)]
        sum = format["add"](terms)
        if sum: sum = format["grouping"](sum)
        if sum: f1 = [sum]
        else: f1 = []
    
        # Compute product of all factors
        return format["multiply"]([f for f in [format["scale factor"]] + f0 + f1])
