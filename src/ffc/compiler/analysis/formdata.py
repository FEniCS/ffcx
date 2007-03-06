__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2005-03-15 -- 2007-03-05"
__copyright__ = "Copyright (C) 2005-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# Modified by Garth N. Wells 2006

# Python modules
import numpy
import sets

# FFC common modules
from ffc.common.debug import *

# FFC fem modules
from ffc.fem.dofmap import *

# FFC language modules
from ffc.compiler.language.integral import *

class FormData:
    """This class holds meta data for a form. The following attributes
    are extracted and stored for a given form:

        form                         - the form
        name                         - the name of the form
        signature                    - the signature of the form
        rank                         - the rank (arity) of the form
        num_coefficients             - the number of coefficients
        num_arguments                - the sum of rank and num_coefficients
        num_terms                    - the number of terms
        num_cell_integrals           - the number of cell integrals
        num_exterior_facet_integrals - the number of exterior facet integrals
        num_interior_facet_integrals - the number of interior facet integrals
        elements                     - the finite elements associated with the form
        dof_maps                     - the dof maps associated with the form
        cell_dimension               - the dimension of the cell

    It is assumed that the indices of the given form have been reassigned.
    """

    def __init__(self, form, name):
        "Create form data for form"

        debug("Extracting form data...")

        self.form                         = form
        self.name                         = name
        self.signature                    = self.__extract_signature(form)
        self.rank                         = self.__extract_rank(form)
        self.num_coefficients             = self.__extract_num_coefficients(form)
        self.num_arguments                = self.rank + self.num_coefficients
        self.num_terms                    = self.__extract_num_terms(form)
        self.num_cell_integrals           = self.__extract_num_cell_integrals(form)
        self.num_exterior_facet_integrals = self.__extract_num_exterior_facet_integrals(form)
        self.num_interior_facet_integrals = self.__extract_num_interior_facet_integrals(form)
        self.elements                     = self.__extract_elements(form, self.rank, self.num_coefficients)
        self.dof_maps                     = self.__extract_dof_maps(self.elements)
        self.cell_dimension               = self.__extract_cell_dimension(self.elements)

        debug("done")

    def __extract_signature(self, form):
        "Extract the signature"
        return str(form)

    def __extract_rank(self, form):
        "Extract the rank"
        return max_index(form, Index.PRIMARY) + 1

    def __extract_num_coefficients(self, form):
        "Extract the number of coefficients"
        return max_index(form, Index.FUNCTION) + 1

    def __extract_num_terms(self, form):
        "Extract the number of terms"
        return len(form.monomials)

    def __extract_num_cell_integrals(self, form):
        "Extract the number of cell integrals"
        integrals = [monomial.integral for monomial in form.monomials if monomial.integral.type == Integral.CELL]
        return self.__extract_num_sub_domains(integrals)

    def __extract_num_exterior_facet_integrals(self, form):
        "Extract the number of exterior facet integrals"
        integrals = [monomial.integral for monomial in form.monomials if monomial.integral.type == Integral.EXTERIOR_FACET]
        return self.__extract_num_sub_domains(integrals)

    def __extract_num_interior_facet_integrals(self, form):
        "Extract the number of interiof facet integrals"
        integrals = [monomial.integral for monomial in form.monomials if monomial.integral.type == Integral.INTERIOR_FACET]
        return self.__extract_num_sub_domains(integrals)
    
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

    def __extract_cell_dimension(self, elements):
        "Extract cell dimension"

        # Should be the same so pick first
        return elements[0].cell_dimension()

    def __extract_num_sub_domains(self, integrals):
        "Extract number of sub domains from list of integrals"
        sub_domains = sets.Set([integral.sub_domain for integral in integrals])
        if len(sub_domains) == 0:
            return 0
        if not (max(sub_domains) + 1 == len(sub_domains) and min(sub_domains) == 0):
            raise FormError, "Sub domains must be numbered from 0 to n - 1"
        return len(sub_domains)

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
