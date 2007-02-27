__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-02-05 -- 2007-02-05"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC common modules
from ffc.common.debug import *

# FFC language modules
from ffc.compiler.language.integral import *

# FFC tensor representation modules
from factorization import *
from referencetensor import *
from geometrytensor import *

class Term:
    """This class represents a tensor contraction A0 : (G0 + G1 + ...)
    of a reference tensor A0 and a sum of geometry tensors G0, G1, ..."""

    def __init__(self, A0, G):
        "Create term A0 : (G0 + G1 + ...)"
        self.A0 = A0
        self.G  = G

class TensorRepresentation:
    """This class represents a given multilinear form as a tensor
    contraction, or more precisely, a sum of tensor contractions for
    each type of integral: cell, exterior facet and interior facet."""

    def __init__(self, form_data):
        "Create tensor representation for given form"

        # Extract form
        form = form_data.form

        # Compute representation of cell tensor
        self.cell_tensor = self.__compute_cell_tensor(form)
        
        # Compute representation of exterior facet tensors
        self.exterior_facet_tensors = self.__compute_exterior_facet_tensors(form)

        # Compute representation of interior facet tensors
        self.interior_facet_tensors = self.__compute_interior_facet_tensors(form)
        
    def __compute_cell_tensor(self, form):
        "Compute representation of cell tensor"
        debug_begin("Computing cell tensor")

        # Check if there are any terms to compute
        if not self.__count_terms(form, Integral.CELL) > 0:
            debug_end()
            return []

        # Compute sum of tensor representations
        terms = self.__compute_terms(form, Integral.CELL, None, None)
        
        debug_end()
        return terms

    def __compute_exterior_facet_tensors(self, form):
        "Compute representation of exterior facet tensors"
        debug_begin("Computing exterior facet tensors")

        # Check if there are any terms to compute
        if not self.__count_terms(form, Integral.EXTERIOR_FACET) > 0:
            debug_end()
            return []

        # Get the number of facets
        num_facets = form.monomials[0].basisfunctions[0].element.num_facets()
        
        # Compute sum of tensor representations for each facet
        terms = [None for i in range(num_facets)]
        for i in range(num_facets):
            terms[i] = self.__compute_terms(form, Integral.EXTERIOR_FACET, i, None)

        debug_end()
        return terms

    def __compute_interior_facet_tensors(self, form):
        "Compute representation of interior facet tensors"
        debug_begin("Computing interior facet tensors")

        # Check if there are any terms to compute
        if not self.__count_terms(form, Integral.INTERIOR_FACET) > 0:
            debug_end()
            return []

        # Get the number of facets
        num_facets = form.monomials[0].basisfunctions[0].element.num_facets()
        
        # Compute sum of tensor representations for each facet-facet combination
        terms = [[None for j in range(num_facets)] for i in range(num_facets)]
        for i in range(num_facets):
            for j in range(num_facets):
                terms[i][j] = self.__compute_terms(form, Integral.INTERIOR_FACET, i, j)

        debug_end()
        return terms

    def __count_terms(self, form, integral_type):
        "Count the number of terms of given type"
        
        num_terms = sum([1 for m in form.monomials if m.integral.type == integral_type])
        debug("Number of terms to consider: %d" % num_terms)
        return num_terms

        # Return empty list if there are no terms to compute
        if num_terms == 0:
            return []

    def __compute_terms(self, form, integral_type, facet0, facet1):
        "Compute terms and factorize common reference tensors"

        # Compute factorization
        debug("Computing factorization...")
        factorization = factorize(form)
        debug("done")

        # Check how many terms we need to compute
        num_terms = sum([1 for m in factorization if m == None])
        debug("Number of terms to compute: %d" % num_terms)

        # Compute terms
        terms = [None for i in range(len(form.monomials))]
        for i in range(len(form.monomials)):

            # Get monomial
            m = form.monomials[i]

            # Only consider monomials of given integral type
            if not m.integral.type == integral_type:
                continue
            
            # Compute geometry tensor for current monomial
            G = GeometryTensor(m)
            
            # Compute reference tensor if not factorized
            if factorization[i] == None:
                # Compute new reference tensor
                A0 = ReferenceTensor(m, facet0, facet1)
                terms[i] = Term(A0, [G])
            else:
                # Add geometry tensor to previous term
                terms[factorization[i]].G += [G]

        # Remove terms not computed (factorized)
        [terms.remove(None) for i in range(len(terms)) if None in terms]

        return terms
