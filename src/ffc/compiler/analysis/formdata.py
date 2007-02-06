__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2005-03-15 -- 2007-02-06"
__copyright__ = "Copyright (C) 2005-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# Modified by Garth N. Wells 2006

# Python modules
import numpy

# FFC common modules
from ffc.common.debug import *

# FFC fem modules
from ffc.fem.dofmap import *

class FormData:
    """This class holds meta data for a form. The following attributes
    are extracted and stored for a given form:

        form             - the form
        rank             - the rank (arity) of the form
        num_coefficients - the number of coefficients
        elements         - the finite elements associated with the form
        dof_maps         - the dof maps associated with the form

    It is assumed that the indices of the given form have been reassigned.
    """

    def __init__(self, form):
        "Create form data for form"

        debug("Extracting form data...")

        self.form = form
        self.rank = self.__extract_rank(form)
        self.num_coefficients = self.__extract_num_coefficients(form)
        self.elements = self.__extract_elements(form, self.rank, self.num_coefficients)
        self.dof_maps = self.__extract_dof_maps(self.elements)

        debug("done")

    def __extract_rank(self, form):
        "Extract the rank"
        return max_index(form, Index.PRIMARY) + 1

    def __extract_num_coefficients(self, form):
        "Extract the number of coefficients"
        return max_index(form, Index.FUNCTION) + 1
    
    def __extract_elements(self, form, rank, num_coefficients):
        """Extract all elements associated with form. The list of elements is
        ordered first by the indices of the corresponding basis functions
        and then by the indices of the corresponding functions."""

        elements = []
        
        # Only need to check first monomial for primary indices
        monomial = form.monomials[0]
        for i in range(rank):
            for v in monomial.basisfunctions:
                if v.index.type == Index.PRIMARY and v.index.index == i:
                    elements += [v.element]

        # Now extract elements corresponding to coefficients
        for i in range(num_coefficients):
            for m in form.monomials:
                for c in m.coefficients:
                    if c.n0.index == i and len(elements) < rank + i + 1:
                        elements += [c.e0]

        # Check that we found an element for each function
        if not len(elements) == rank + num_coefficients:
            raise FormError, (form, "Unable to extract all elements.")

        return elements

    def __extract_dof_maps(self, elements):
        "Extract (generate) dof maps for all elements"
        return [DofMap(element) for element in elements]

    def __repr__(self):
        "Pretty print"
        return """\
Form:         %s
Rank:         %s
Coefficients: %s
Elements:     %s
Dof maps:     %s""" % (str(self.form),
                       str(self.rank),
                       str(self.num_coefficients),
                       str(self.elements),
                       str(self.dof_maps))
