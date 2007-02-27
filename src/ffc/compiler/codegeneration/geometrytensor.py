"Code generation for geometry tensor (for tensor representation)"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-03 -- 2007-02-27"
__copyright__ = "Copyright (C) 2004-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC language modules
from ffc.compiler.language.index import *

def generate_geometry_tensor(terms, format):
    "Generate list of declarations for computation of geometry tensors"

    # Generate code as a list of declarations
    code = []    

    # Iterate over all terms
    for j in range(len(terms)):
        
        # Get list of secondary indices (should be the same so pick first)
        aindices = terms[j].G[0].a.indices
        if not aindices: aindices = [[]]

        # Iterate over secondary indices
        for a in aindices:
            
            # Sum factorized values
            name = format["geometry tensor"](j, a)
            value = format["sum"]([__generate_entry(G, a, format) for G in terms[j].G])

            # Add declaration
            code += [(name, value)]

    return code

def __generate_entry(G, a, format):
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
            factors += [format["transform"](t.index0([], a, [], []), \
                                            t.index1([], a, [], []), \
                                            t.restriction),]
    monomial = format["multiplication"](factors)
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
                factors += [format["transform"](t.index0([], a, [], b), \
                                                t.index1([], a, [], b), \
                                                t.restriction)]
        terms += [format["multiplication"](factors)]
    sum = format["sum"](terms)
    if sum: sum = format["grouping"](sum)
    if sum: f1 = [sum]
    else: f1 = []

    # Compute product of all factors
    return format["multiplication"]([f for f in [format["determinant"]] + f0 + f1])
