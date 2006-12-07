__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2006-12-01 -- 2006-12-04"
__copyright__ = "Copyright (C) 2006 Anders Logg"
__license__  = "GNU GPL Version 2"

# Python modules
import numpy

# FFC common modules
from ffc.common.debug import *
from ffc.common.constants import *
from ffc.common.exceptions import *

# FFC compiler modules
from tensorrepresentation import *
from multiindex import *
from integral import *

class InteriorFacetTensor:
    """An InteriorElementTensor represents the interior facet tensor for
    a multi-linear form and consists of a list of Terms, each containing
    a pair of a ReferenceTensor and a GeometryTensor.

    Attributes:

        terms   - a list of Terms (products A0 * GS)
        aS      - a list of precomputed element tensor declarations
        a0      - a list of precomputed reference tensor declarations
        gS      - a list of precomputed geometry tensor declarations
        facet0  - local number of facet 0 ('+' facet for interior facet tensor)
        facet1  - local number of facet 1 ('-' facet for interior facet tensor)
        num_ops - number of operations in computation of element tensor
    """

    def __init__(self, sum, format, cS_used, gS_used, options, facet0, facet1, alignment):
        "Create ElementTensor."

        print "Compiling element tensor: " + str(sum) + " " + str((facet0, facet1, alignment))

        # Reset number of operations
        self.num_ops = 0

        # Compute terms
        self.terms = compute_terms(sum, Integral.INTERIOR_FACET, facet0, facet1, alignment)
        if len(self.terms) == 0:
            return

        # Reorder entries to compute reference tensor from reduced reference tensor
        self.__reorder_entries()

        # Compute element tensor declarations
        (self.aS, self.num_ops) = compute_element_tensor(self.terms, format, options)

        # Compute reference tensor declarations
        self.a0 = compute_reference_tensor(self.terms, format)
        
        # Compute geometry tensor declarations
        check_used(self.terms, format, gS_used)
        self.gS = compute_geometry_tensor(self.terms, format, gS_used, cS_used)

        # Save facets
        self.facet0 = facet0
        self.facet1 = facet1

    def __reorder_entries(self):
        """Reorder entries in the reference tensor to compute the reference tensor
        from the reduced reference tensor."""

        for term in self.terms:

            # Compute restrictions corresponding to indices
            restrictions = self.__compute_restrictions(term)

            # Compute size of reference tensor
            dims = term.A0.i.dims + term.A0.a.dims
            new_idims = [2*d for d in term.A0.i.dims]
            new_adims = [2*d for d in term.A0.a.dims]
            new_dims = new_idims + new_adims

            # Compute position where to insert
            position = []
            for i in range(len(restrictions)):
                dim = dims[i]
                if restrictions[i] == Restriction.PLUS:
                    position = position + [slice(0, dim)]
                elif restrictions[i] == Restriction.MINUS:
                    position = position + [slice(dim, 2*dim)]
                else:
                    raise FormError, (term.product, "Missing restriction for function in facet integral.")

            # Initialize empty reference tensor of double size in each dimension
            A0 = numpy.zeros(new_dims, dtype= numpy.float)

            # Insert reduced reference tensor into reference tensor
            A0[position] = term.A0.A0
            term.A0.A0 = A0

            print A0

            # Reinitialize indices to new size
            term.A0.i = MultiIndex(new_idims)
            term.A0.a = MultiIndex(new_adims)

    def __compute_restrictions(self, term):
        """Compute restrictions corresponding to indices for given term."""

        # Get dimensions for primary and secondary indices
        idims = term.A0.i.dims
        adims = term.A0.a.dims
        
        # Get basis functions for term
        basisfunctions = term.product.basisfunctions
        
        # Create empty list of restrictions for indices
        restrictions = [None for i in range(len(idims) + len(adims))]
        
        # Extract restrictions corresponding to primary indices
        for i in range(len(idims)):
            for v in basisfunctions:
                if v.index.type == Index.PRIMARY and v.index.index == i:
                    restrictions[i] = v.restriction
                    break

        # Extract restrictions corresponding to primary indices
        for i in range(len(adims)):
            for v in basisfunctions:
                if v.index.type == Index.SECONDARY and v.index.index == i:
                    restrictions[len(idims) + i] = v.restriction
                    break

        # Check that we got the restrictions
        for restriction in restrictions:
            if restriction == None:
                raise FormError, (term.product, "Missing restriction for function in facet integral.")

        return restrictions
