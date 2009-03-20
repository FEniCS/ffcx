__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-02-05 -- 2009-03-05"
__copyright__ = "Copyright (C) 2007-2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# UFL modules
from ufl.classes import Form
from ufl.integral import Measure
from ufl.algorithms import extract_basis_functions

# FFC common modules
from ffc.common.log import debug, info, warning, error, begin, end, set_level, INFO

# FFC fem modules
from ffc.fem import create_element 

# FFC language modules
#from ffc.compiler.language.integral import *

# FFC tensor representation modules
from monomialextraction import extract_monomial_form, MonomialForm
from monomialtransformation import transform_monomial_form
from uflreferencetensor import ReferenceTensor
from uflgeometrytensor import GeometryTensor
from ufltensorreordering import reorder_entries

class TensorContraction:
    "This class represents a tensor contraction A^K = A^0 : G_K."

    def __init__(self, monomial, domain_type, facet0=None, facet1=None):
        "Create tensor contraction for given monomial."
        self.monomial = monomial
        self.A0 = ReferenceTensor(monomial, domain_type, facet0, facet1)
        # FIXME: Does not need to be list if we find that factorization is not needed
        self.GK = [GeometryTensor(monomial)]

class TensorRepresentation:
    """This class represents a given multilinear form as a tensor
    contraction, or more precisely, a sum of tensor contractions for
    each type of integral: cell, exterior facet and interior facet.

    Attributes:

        cell_integrals           - list of list of terms for sub domains
        exterior_facet_integrals - list of list of list of terms for sub domains and facets
        interior_facet_integrals - list of list of list of list of terms for sub domains and facet combinations
        geometric_dimension      - geometric dimension of form
        num_facets               - number of cell facets

    Each term is represented as a TensorContraction.
    """

    def __init__(self, form_data):
        "Create tensor representation for given form."

        print "Input form:", form_data.form

        # Extract integrals integrals for tensor representation
        form = _extract_tensor_integrals(form_data.form)

        print "Extracted form:", form

        # FIXME: Temporary fix
        if len(form.integrals()) == 0:
            self.cell_integrals = []
            self.interior_facet_integrals = []
            self.exterior_facet_integrals = []
            return

        # Extract monomial representation
        monomial_form = extract_monomial_form(form)
        print ""
        print "Monomial form:", monomial_form

        # Transform monomial form to reference element
        transform_monomial_form(monomial_form)
        print "Transformed monomial:", monomial_form

        # Compute representation of cell tensor
        n = form_data.num_cell_integrals
        self.cell_integrals = [_compute_cell_tensor(monomial_form, form_data, i) for i in range(n)]
        
        # Compute representation of exterior facet tensors
        n = form_data.num_exterior_facet_integrals
        self.exterior_facet_integrals = [_compute_exterior_facet_tensors(monomial_form, form_data, i) for i in range(n)]

        # Compute representation of interior facet tensors
        n = form_data.num_interior_facet_integrals
        self.interior_facet_integrals = [_compute_interior_facet_tensors(monomial_form, form_data, i) for i in range(n)]

        # Extract form data needed by code generation
        self.geometric_dimension = form_data.geometric_dimension
        self.num_facets = form_data.num_facets
        
def _extract_tensor_integrals(form):
    "Extract form containing only tensor representation integrals."

    new_form = Form([])
    for integral in form.integrals():
        if integral.measure().metadata()["ffc_representation"] == "tensor":
            new_form += Form([integral])
    return new_form

def _compute_cell_tensor(monomial_form, form_data, sub_domain):
    "Compute representation of cell tensor."
    
    begin("Computing cell tensor")
    
    # Extract cell integrals
    monomial_form = _extract_integrals(monomial_form, form_data, Measure.CELL, sub_domain)
    
    # Compute sum of tensor representations
    terms = _compute_terms(monomial_form, Measure.CELL, None, None)

    end()
    
    return terms

def _compute_exterior_facet_tensors(monomial_form, form_data, sub_domain):
    "Compute representation of exterior facet tensors."

    begin("Computing exterior facet tensors")

    # Extract exterior facet integrals
    monomial_form = _extract_integrals(monomial_form, form_data, Measure.EXTERIOR_FACET, sub_domain)
    
    # Compute sum of tensor representations for each facet
    terms = [None for i in range(form_data.num_facets)]
    for i in range(form_data.num_facets):
        terms[i] = _compute_terms(monomial_form, Measure.EXTERIOR_FACET, i, None)

    end()

    return terms

def _compute_interior_facet_tensors(monomial_form, form_data, sub_domain):
    "Compute representation of interior facet tensors."

    begin("Computing interior facet tensors")

    # Extract interior facet integrals
    monomial_form = _extract_integrals(monomial_form, form_data, Measure.INTERIOR_FACET, sub_domain)
    
    # Compute sum of tensor representations for each facet-facet combination
    terms = [[None for j in range(form_data.num_facets)] for i in range(form_data.num_facets)]
    for i in range(form_data.num_facets):
        for j in range(form_data.num_facets):
            terms[i][j] = _compute_terms(monomial_form, Measure.INTERIOR_FACET, i, j)
            reorder_entries(terms[i][j])

    end()

    return terms

def _extract_integrals(monomial_form, form_data, domain_type, sub_domain):
    "Extract subset of form matching given domain type."

    new_form = MonomialForm()
    for (integrand, measure) in monomial_form:
        if measure.domain_type() == domain_type and measure.domain_id() == sub_domain:
            new_form.append(integrand, measure)
    return new_form

def _compute_terms(monomial_form, domain_type, facet0, facet1):
    "Compute list of tensor contraction terms for monomial form."

    # Compute terms
    terms = []
    for (integrand, measure) in monomial_form:

        # Only consider monomials of given integral type
        if not measure.domain_type() == domain_type:
            continue

        # Iterate over monomials of integrand
        for monomial in integrand.monomials:
            terms.append(TensorContraction(monomial, domain_type, facet0, facet1))

    return terms
