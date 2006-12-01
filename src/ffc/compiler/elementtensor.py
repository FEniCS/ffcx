__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-06 -- 2006-12-01"
__copyright__ = "Copyright (C) 2004-2006 Anders Logg"
__license__  = "GNU GPL Version 2"

# Modified by Garth N. Wells 2006

# FFC common modules
from ffc.common.debug import *
from ffc.common.constants import *
from ffc.common.exceptions import *

# FFC compiler modules
from term import *
from reorder import *
from declaration import *
from tensorrepresentation import *

class ElementTensor:
    """An ElementTensor represents the element tensor of a
    multi-linear form and consists of a list of Terms, each containing
    a pair of a ReferenceTensor and a GeometryTensor.

    Attributes:

        terms   - a list of Terms (products A0 * GK)
        a0      - a list of precomputed reference tensor declarations
        gK      - a list of precomputed geometry tensor declarations
        aK      - a list of precomputed element tensor declarations
        facet0  - local number of facet 0 ('+' facet for interior facet tensor)
        facet1  - local number of facet 1 ('-' facet for interior facet tensor)
        num_ops - number of operations in computation of element tensor

    Note that this class is used both for the representation of
    standard element tensor A^K and exterior/interior facet tensors A^S.
    """

    def __init__(self, sum, type, format, cK_used, gK_used, options, facet0, facet1):
        "Create ElementTensor."

        # Reset number of operations
        self.num_ops = 0

        # Check if there are any terms to compute
        num_terms = terms_to_compile(sum, type)
        debug("Number of terms to compile: %d" % num_terms)
        if num_terms == 0:
            self.terms = []
            return

        # Reorder indices and compute factorization
        factorization = reorder_indices(sum)

        # Compute terms
        self.terms = [None for i in range(len(sum.products))]
        for i in range(len(sum.products)):
            debug("Compiling term %d" % i, 1)
            p = sum.products[i]
            if p.integral.type == type:
                # Compute geometry tensor
                GK = GeometryTensor(p)
                # Check if reference tensor should be computed
                if factorization[i] == None:
                    # Compute reference tensor and add term
                    A0 = ReferenceTensor(p, facet0, facet1)
                    self.terms[i] = Term(p, A0, [GK])
                else:
                    # Add geometry tensor to previous term
                    self.terms[factorization[i]].G += [GK]
        debug("All terms compiled", 1)

        # Remove terms not computed (factorized)
        [self.terms.remove(None) for i in range(len(self.terms)) if None in self.terms]
    
        # Compute reference tensor declarations (just need to pick the values)
        self.a0 = compute_reference_tensor(self.terms, format)

        # Compute element tensor declarations
        self.aK = compute_element_tensor(self.terms, format, options)

        # Compute geometry tensor declarations
        check_used(self.terms, format, gK_used)
        self.gK = compute_geometry_tensor(self.terms, format, gK_used, cK_used)

        # Save facets
        self.facet0 = facet0
        self.facet1 = facet1

        return
