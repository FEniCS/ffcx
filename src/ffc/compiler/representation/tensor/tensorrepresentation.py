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

    def __init__(self, A0, Gs):
        "Create term A0 : (G0 + G1 + ...)"
        self.A0 = A0
        self.Gs = Gs

class TensorRepresentation:
    """This class represents a given multilinear form as a tensor
    contraction, or more precisely, a sum of tensor contractions for
    each type of integral: cell, exterior facet and interior facet."""

    def __init__(self, form_data):
        "Create tensor representation for given form"

        # Extract form
        form = form_data.form

        # Compute representation of cell tensor
        self.__cell_tensor = self.__compute_cell_tensor(form)
        
        # Compute representation of exterior facet tensors
        self.__exterior_facet_tensors = self.__compute_exterior_facet_tensors(form)

        # Compute representation of interior facet tensors
        self.__interior_facet_tensors = self.__compute_interior_facet_tensors(form)
        
    def cell_tensor():
        "Return representation of cell tensor"
        return self.__cell_tensor

    def exterior_facet_tensor(self, i, j):
        "Return representation of exterior facet tensor on local facet i"
        return self.__exterior_facet_tensors[i]

    def interior_facet_tensor(self, i):
        """Return representation of interior facet tensor on intersecton
        of local facets i and j"""
        
        return self.__interior_facet_tensors[i][j]

    def __compute_cell_tensor(self, form):
        "Compute representation of cell tensor"
        debug_begin("Computing cell tensor")
        
        terms = self.__compute_terms(form, Integral.CELL, None, None)

        debug_end()

        return "Not implemented"

    def __compute_exterior_facet_tensors(self, form):
        "Compute representation of exterior facet tensors"
        debug_begin("Computing exterior facet tensors")

        debug("Not implemented")

        debug_end()
        
        return "Not implemented"

    def __compute_interior_facet_tensors(self, form):
        "Compute representation of interior facet tensors"
        debug_begin("Computing interior facet tensors")

        debug("Not implemented")
        
        debug_end()
        
        return "Not implemented"

    def __compute_terms(self, form, integral_type, facet0, facet1):
        "Compute terms and factorize common reference tensors"

        # Count the number of terms to compute
        num_terms = sum([1 for m in form.monomials if m.integral.type == integral_type])
        debug("Number of terms to consider: %d" % num_terms)

        # Return empty list if there are no terms to compute
        if num_terms == 0:
            return []

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
                A0 = ReferenceTensor(m, facet0, facet1)
                terms[i] = Term(A0, [G])
            else:
                # Add geometry tensor to previous term
                terms[factorization[i]].Gs += [G]

        # Remove terms not computed (factorized)
        [terms.remove(None) for i in range(len(terms)) if None in terms]
