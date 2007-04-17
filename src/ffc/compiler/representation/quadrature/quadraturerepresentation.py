__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2007-03-16 -- 2007-03-23"
__copyright__ = "Copyright (C) 2007 Kristian B. Oelgaard"
__license__  = "GNU GPL Version 2"

# FFC common modules
from ffc.common.debug import *

# FFC language modules
from ffc.compiler.language.integral import *

# FFC quadrature representation modules
from elementtensor import *

# obsolete modules?
#from monomialintegration import *

#from factorization import *
#from referencetensor import *
#from geometrytensor import *
#from tensorreordering import *

class QuadratureRepresentation:
    """This class uses quadrature to represent a given multilinear form.

    Attributes:

        form                   - the form generating the quadrature representation
        cell_tensor            - the representation of the cell tensor
        exterior_facet_tensors - the representation of the interior facet tensors,
                                 one for each facet
        interior_facet_tensors - the representation of the exterior facet tensors,
                                 one for each facet-facet combination
    """

    def __init__(self, form_data):
        "Create tensor representation for given form"

        # Extract form
        form = form_data.form

#        print "\n quadraturerepresentations __init__, form: \n", form

        # Save form
        self.form = form

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
#        print "monomials[0]: ", monomials[0]
#        print "monomials[1]: ", monomials[1]

        # Compute factorization
        #FIXME: this will mean something else for quadrature
        factorization = self.__compute_factorization(monomials)

        # Compute sum of tensor representations
        tensors = self.__compute_tensors(monomials, factorization, Integral.CELL, None, None)

        debug_end()
        return tensors

    # FIXME: not implemented/tested
    def __compute_exterior_facet_tensors(self, form):
        "Compute representation of exterior facet tensors"
        debug_begin("Computing exterior facet tensors")

        # Extract monomials
        monomials = self.__extract_monomials(form, Integral.EXTERIOR_FACET)
        if len(monomials) == 0:
            debug_end()
            return []

        # Compute factorization
#        factorization = self.__compute_factorization(monomials)

        # Get the number of facets
#        num_facets = form.monomials[0].basisfunctions[0].element.num_facets()

#        debug("Number of facets to consider: %d" % num_facets)
        
        # Compute sum of tensor representations for each facet
#        terms = [None for i in range(num_facets)]
#        for i in range(num_facets):
#            terms[i] = self.__compute_tensors(monomials, factorization, Integral.EXTERIOR_FACET, i, None)

        tensors = []
        debug_end()
        return tensors

    # FIXME: not implemented/tested
    def __compute_interior_facet_tensors(self, form):
        "Compute representation of interior facet tensors"
        debug_begin("Computing interior facet tensors")

        # Extract monomials
        monomials = self.__extract_monomials(form, Integral.INTERIOR_FACET)
        if len(monomials) == 0:
            debug_end()
            return []

        # Compute factorization
#        factorization = self.__compute_factorization(monomials)

        # Get the number of facets
#        num_facets = form.monomials[0].basisfunctions[0].element.num_facets()

#        debug("Number of facets to consider: %d x %d" % (num_facets, num_facets))
        
        # Compute sum of tensor representations for each facet-facet combination
#        terms = [[None for j in range(num_facets)] for i in range(num_facets)]
#        for i in range(num_facets):
#            for j in range(num_facets):
#                terms[i][j] = self.__compute_tensors(monomials, factorization, Integral.INTERIOR_FACET, i, j)
#                reorder_entries(terms[i][j])

        tensors = []               
        debug_end()
        return tensors

    def __extract_monomials(self, form, integral_type):
        "Extract monomials"

        # Extract monomials of given type
        monomials = [m for m in form.monomials if m.integral.type == integral_type]
        if len(monomials) > 0:
            debug("Number of terms to consider: %d" % len(monomials))
        else:
            debug("No terms")

        return monomials

    #FIXME: this will mean something else for quadrature, returns None
    def __compute_factorization(self, monomials):
        "Compute factorization"

#        factorization = factorize(monomials)
#        print "factor"
#        print "len(monomials): ", len(monomials)
        factorization = [None for i in range(len(monomials))]

        num_terms = sum([1 for m in factorization if m == None])
        debug("Number of terms to compute: %d" % num_terms)

        return factorization

    def __compute_tensors(self, monomials, factorization, integral_type, facet0, facet1):
        "Compute terms and factorize common reference tensors"

        # Compute terms
        num_tensors = len(monomials)
        tensors = [None for i in range(num_tensors)]

        for i in range(num_tensors):

            # Get monomial
            m = monomials[i]

#            print "compute terms, m: ", m
#            print "compute terms, m.basisfunctions: ", m.basisfunctions

#            print "degree: ", compute_degree(m.basisfunctions)
#            print "num_gauss_points: ", num

            tensors[i] = ElementTensor(m, facet0, facet1)

            # Only consider monomials of given integral type
#            if not m.integral.type == integral_type:
#                continue
            
            # Compute geometry tensor for current monomial
#            G = GeometryTensor(m)
            
            # Compute reference tensor if not factorized
#            self.__debug(i, facet0, facet1)
#            if factorization[i] == None:
                # Compute new reference tensor
#                A0 = ReferenceTensor(m, facet0, facet1)
#                terms[i] = Term(m, A0, [G])
#                debug("done")
#            else:
                # Add geometry tensor to previous term
#                terms[factorization[i]].G += [G]
#                debug("factorized")

        # Remove terms not computed (factorized)
#        [terms.remove(None) for i in range(len(terms)) if None in terms]

        return tensors

    def __debug(self, i, facet0, facet1):
        "Fancy printing of progress"
        if facet0 == facet1 == None:
            debug("Computing quadrature representation for term %d..." % i)
        elif facet1 == None:
            debug("Computing quadrature representation for facet %d, term %d..." % (facet0, i))
        else:
            debug("Computing quadrature representation for facets (%d, %d), term %d..." % (facet0, facet1, i))
