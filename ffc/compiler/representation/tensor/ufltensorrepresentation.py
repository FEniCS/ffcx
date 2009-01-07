__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-02-05 -- 2007-03-09"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# UFL modules
from ufl.algorithms import extract_monomials

# FFC common modules
from ffc.common.debug import *

# FFC language modules
from ffc.compiler.language.integral import *

# FFC tensor representation modules
from factorization import *
from referencetensor import *
from geometrytensor import *
from tensorreordering import *

class Term:
    """This class represents a tensor contraction A0 : (G0 + G1 + ...)
    of a reference tensor A0 and a sum of geometry tensors G0, G1, ..."""

    def __init__(self, monomial, A0, G):
        "Create term A0 : (G0 + G1 + ...)"
        self.monomial = monomial
        self.A0 = A0
        self.G  = G

class TensorRepresentation:
    """This class represents a given multilinear form as a tensor
    contraction, or more precisely, a sum of tensor contractions for
    each type of integral: cell, exterior facet and interior facet.

    Attributes:

        form                   - the form generating the tensor contraction
        cell_tensor            - the representation of the cell tensor
        exterior_facet_tensors - the representation of the interior facet tensors,
                                 one for each facet
        interior_facet_tensors - the representation of the exterior facet tensors,
                                 one for each facet-facet combination
    """

    def __init__(self, form_data, num_quadrature_points):
        "Create tensor representation for given form"

        print "in constructor"

        # Extract form
        form = form_data.form

        print form

        # Convert to monomial representation
        monomials = extract_monomials(form)

        print monomials

        return

        # Compute representation of cell tensor
        self.cell_tensor = self.__compute_cell_tensor(form)
        
        # Compute representation of exterior facet tensors
        self.exterior_facet_tensors = self.__compute_exterior_facet_tensors(form)

        # Compute representation of interior facet tensors
        self.interior_facet_tensors = self.__compute_interior_facet_tensors(form)
        
    def __compute_cell_tensor(self, form):
        "Compute representation of cell tensor"
        debug_begin("Computing cell tensor")

        # Extract monomials
        monomials = self.__extract_monomials(form, Integral.CELL)
        if len(monomials) == 0:
            debug_end()
            return []

        # Compute factorization
        factorization = self.__compute_factorization(monomials)

        # Compute sum of tensor representations
        terms = self.__compute_terms(monomials, factorization, Integral.CELL, None, None)
        
        debug_end()
        return terms

    def __compute_exterior_facet_tensors(self, form):
        "Compute representation of exterior facet tensors"
        debug_begin("Computing exterior facet tensors")

        # Extract monomials
        monomials = self.__extract_monomials(form, Integral.EXTERIOR_FACET)
        if len(monomials) == 0:
            debug_end()
            return []

        # Compute factorization
        factorization = self.__compute_factorization(monomials)

        # Get the number of facets
        num_facets = form.monomials[0].basisfunctions[0].element.num_facets()

        debug("Number of facets to consider: %d" % num_facets)
        
        # Compute sum of tensor representations for each facet
        terms = [None for i in range(num_facets)]
        for i in range(num_facets):
            terms[i] = self.__compute_terms(monomials, factorization, Integral.EXTERIOR_FACET, i, None)

        debug_end()
        return terms

    def __compute_interior_facet_tensors(self, form):
        "Compute representation of interior facet tensors"
        debug_begin("Computing interior facet tensors")

        # Extract monomials
        monomials = self.__extract_monomials(form, Integral.INTERIOR_FACET)
        if len(monomials) == 0:
            debug_end()
            return []

        # Compute factorization
        factorization = self.__compute_factorization(monomials)

        # Get the number of facets
        num_facets = form.monomials[0].basisfunctions[0].element.num_facets()

        debug("Number of facets to consider: %d x %d" % (num_facets, num_facets))
        
        # Compute sum of tensor representations for each facet-facet combination
        terms = [[None for j in range(num_facets)] for i in range(num_facets)]
        for i in range(num_facets):
            for j in range(num_facets):
                terms[i][j] = self.__compute_terms(monomials, factorization, Integral.INTERIOR_FACET, i, j)
                reorder_entries(terms[i][j])
                
        debug_end()
        return terms

    def __extract_monomials(self, form, integral_type):
        "Extract monomials and factorize"

        # Extract monomials of given type
        monomials = [m for m in form.monomials if m.integral.type == integral_type]
        if len(monomials) > 0:
            debug("Number of terms to consider: %d" % len(monomials))
        else:
            debug("No terms")

        return monomials

    def __compute_factorization(self, monomials):
        "Compute factorization"

        factorization = factorize(monomials)
        num_terms = sum([1 for m in factorization if m == None])
        debug("Number of terms to compute: %d" % num_terms)
        return factorization

    def __compute_terms(self, monomials, factorization, integral_type, facet0, facet1):
        "Compute terms and factorize common reference tensors"

        # Compute terms
        terms = [None for i in range(len(monomials))]
        for i in range(len(monomials)):

            # Get monomial
            m = monomials[i]

            # Only consider monomials of given integral type
            if not m.integral.type == integral_type:
                continue
            
            # Compute geometry tensor for current monomial
            G = GeometryTensor(m)
            
            # Compute reference tensor if not factorized
            self.__debug(i, facet0, facet1)
            if factorization[i] == None:
                # Compute new reference tensor
                A0 = ReferenceTensor(m, facet0, facet1)
                terms[i] = Term(m, A0, [G])
                debug("done")
            else:
                # Add geometry tensor to previous term
                terms[factorization[i]].G += [G]
                debug("factorized")

        # Remove terms not computed (factorized)
        [terms.remove(None) for i in range(len(terms)) if None in terms]

        return terms

    def __debug(self, i, facet0, facet1):
        "Fancy printing of progress"
        if facet0 == facet1 == None:
            debug("Computing tensor representation for term %d..." % i)
        elif facet1 == None:
            debug("Computing tensor representation for facet %d, term %d..." % (facet0, i))
        else:
            debug("Computing tensor representation for facets (%d, %d), term %d..." % (facet0, facet1, i))
