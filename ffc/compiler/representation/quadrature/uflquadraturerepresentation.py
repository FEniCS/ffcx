"Quadrature representation class"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2009-01-07 -- 2009-01-07"
__copyright__ = "Copyright (C) 2009 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# FFC common modules
from ffc.common.debug import *

# FFC language modules
from ffc.compiler.language.integral import *

# FFC quadrature representation modules
from elementtensor import *

#from factorization import *
#from tensorreordering import *

class QuadratureRepresentation:
    """This class uses quadrature to represent a given multilinear form.

    Attributes:

        form                            - the form generating the quadrature representation
        cell_tensor                     - the representation of the cell tensor
        exterior_facet_tensors          - the representation of the interior facet tensors,
                                          one for each facet
        interior_facet_tensors          - the representation of the exterior facet tensors,
                                          one for each facet-facet combination
        num_user_specified_quad_points  - the number of desired quadrature points specified
                                          by the user. Will be used for ALL terms

    """

    def __init__(self, form_data, num_quadrature_points):
        "Create tensor representation for given form"

        # Extract form
        form = form_data.form

        # Save form
        self.form = form

        # Set number of specified quadrature points
        self.num_user_specified_quad_points = num_quadrature_points

        print "QR: form", form

        # Compute representation of cell tensor
#        self.cell_tensor = self.__compute_cell_tensor(form)
        
        # Compute representation of exterior facet tensors
#        self.exterior_facet_tensors = self.__compute_exterior_facet_tensors(form)

        # Compute representation of interior facet tensors
#        self.interior_facet_tensors = self.__compute_interior_facet_tensors(form)
        
    def __compute_cell_tensor(self, form):
        "Compute representation of cell tensor"
        debug_begin("Computing cell tensor")

        # Extract monomials
        monomials = self.__extract_monomials(form, Integral.CELL)
        if len(monomials) == 0:
            debug_end()
            return []

        # Compute factorization
        #FIXME: this will mean something else for quadrature
        factorization = self.__compute_factorization(monomials)

        # Compute sum of tensor representations
        tensors = self.__compute_tensors(monomials, factorization, Integral.CELL, None, None)

        debug_end()
        return tensors

    def __compute_exterior_facet_tensors(self, form):
        "Compute representation of exterior facet tensors"
        debug_begin("Computing exterior facet tensors")

        # Extract monomials
        monomials = self.__extract_monomials(form, Integral.EXTERIOR_FACET)
        if len(monomials) == 0:
            debug_end()
            return []

        # Compute factorization
        #FIXME: this will mean something else for quadrature
        factorization = self.__compute_factorization(monomials)

        # Get the number of facets
        num_facets = form.monomials[0].basisfunctions[0].element.num_facets()

        debug("Number of facets to consider: %d" % num_facets)
        
        # Compute sum of tensor representations for each facet
        tensors = [None for i in range(num_facets)]
        for i in range(num_facets):
            tensors[i] = self.__compute_tensors(monomials, factorization, Integral.EXTERIOR_FACET, i, None)

        debug_end()
        return tensors

    def __compute_interior_facet_tensors(self, form):
        "Compute representation of interior facet tensors"
        debug_begin("Computing interior facet tensors")

        # Extract monomials
        monomials = self.__extract_monomials(form, Integral.INTERIOR_FACET)
        if len(monomials) == 0:
            debug_end()
            return []

        # Compute factorization
        #FIXME: this will mean something else for quadrature
        factorization = self.__compute_factorization(monomials)

        # Get the number of facets
        num_facets = form.monomials[0].basisfunctions[0].element.num_facets()

        debug("Number of facets to consider: %d x %d" % (num_facets, num_facets))
        
        # Compute sum of tensor representations for each facet-facet combination
        tensors = [[None for j in range(num_facets)] for i in range(num_facets)]
        for i in range(num_facets):
            for j in range(num_facets):
                tensors[i][j] = self.__compute_tensors(monomials, factorization, Integral.INTERIOR_FACET, i, j)
#                reorder_entries(terms[i][j])

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

            # Only consider monomials of given integral type
            if not m.integral.type == integral_type:
                continue

            # Compute element tensor
            self.__debug(i, facet0, facet1)
            tensors[i] = ElementTensor(m, facet0, facet1, self.num_user_specified_quad_points)
            debug("done")

        return tensors

    def __debug(self, i, facet0, facet1):
        "Fancy printing of progress"
        if facet0 == facet1 == None:
            debug("Computing quadrature representation for term %d..." % i)
        elif facet1 == None:
            debug("Computing quadrature representation for facet %d, term %d..." % (facet0, i))
        else:
            debug("Computing quadrature representation for facets (%d, %d), term %d..." % (facet0, facet1, i))
