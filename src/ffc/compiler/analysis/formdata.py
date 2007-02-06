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

# FFC language modules
from ffc.compiler.language import *

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
        self.elements = self.__extract_elements(form, self.rank)
        self.dof_maps = self.__extract_dof_maps(self.elements)

        debug("done")

    def __extract_rank(self, form):
        "Extract the rank"
        return max_index(form, Index.PRIMARY) + 1

    def __extract_num_coefficients(self, form):
        "Extract the number of coefficients"
        return max_index(form, Index.FUNCTION) + 1
    
    def __extract_elements(self, form, rank):
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
    
def find_elesdfments(form, nfunctions):
    """Return a list of FiniteElements associated with the (original)
    function spaces of the Functions appearing in the form."""

    # List of elements used for functions
    elements = [None for j in range(nfunctions)]

    # Iterate over all Coefficients in all Monomials
    for p in form.monomials:
        for c in p.coefficients:
            elements[c.n0.index] = c.e0

    # Check that we found an element for each function
    for element in elements:
        if not element:
            raise FormError, (form, "Unable to find element for each function.")

    if elements:
        debug("Finite elements for functions: " + str(elements), 0)
          
    return elements

def find_projections(form, nprojections):
    """Return a list of tuples (n0, n1, e0, e1, P) defining the
    projections of all Functions appearing in the form."""

    # List of projections used for functions
    projections = [None for j in range(nprojections)]

    # Iterate over all Coefficients in all Monomials
    for p in form.monomials:
        for c in p.coefficients:
            projections[c.n1.index] = (c.n0.index, c.n1.index, c.e0, c.e1, c.P)

    # Check that we found an element for each projection
    for projection in projections:
        if not projection:
            raise FormError, (form, "Unable to find a projection for each function.")

    return projections
