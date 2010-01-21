"""This module implements the representation of a multilinear form as
a sum of tensor contractions.

The following possible optimizations are currently not implemented but
might be (re-)implemented in a future version of FFC

  1. Factorization of common reference tensors
  2. FErari optimizations
"""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-02-05"
__copyright__ = "Copyright (C) 2007-2010 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian B. Oelgaard, 2010.
# Last changed: 2010-01-21

# UFL modules
from ufl.classes import Form
from ufl.classes import Measure
from ufl.classes import Integral

# FFC modules
from ffc.log import info

# FFC tensor representation modules
from ffc.tensor.monomialextraction import extract_monomial_form
from ffc.tensor.monomialextraction import MonomialForm
from ffc.tensor.monomialtransformation import transform_monomial_form
from ffc.tensor.referencetensor import ReferenceTensor
from ffc.tensor.geometrytensor import GeometryTensor
from ffc.tensor.tensorreordering import reorder_entries

def compute_integral_ir(form, form_data, form_id, options):
    "Compute intermediate represention of form integrals."

    # Extract integrals that should be computed by tensor representation
    form = _extract_tensor_integrals(form, form_data)

    # Check number of integrals
    num_integrals = len(form.integrals())
    if num_integrals == 0: return []

    info("Computing tensor representation")

    # Extract monomial representation
    monomial_form = extract_monomial_form(form, form_data)

    # Transform monomial form to reference element
    transform_monomial_form(monomial_form)
    m = monomial_form

    # List of representations, one for each integral
    ir = []

    # Compute representation of cell tensors
    for i in range(form_data.num_cell_domains):
        AK = _compute_cell_tensor(m, form_data, i)
        ir.append({"AK": AK,
                   "representation": "tensor",
                   "integral_type": "cell_integral",
                   "form_id": form_id,
                   "sub_domain": i,
                   "geometric_dimension": form_data.geometric_dimension,
                   "num_facets": form_data.num_facets})

    # Compute representation of exterior facet tensors
    for i in range(form_data.num_exterior_facet_domains):
        AK = _compute_exterior_facet_tensors(m, form_data, i)
        ir.append({"AK": AK,
                   "representation": "tensor",
                   "integral_type": "exterior_facet_integral",
                   "form_id": form_id,
                   "sub_domain": i,
                   "geometric_dimension": form_data.geometric_dimension,
                   "num_facets": form_data.num_facets})

    # Compute representation of interior facet tensors
    for i in range(form_data.num_interior_facet_domains):
        AK = _compute_interior_facet_tensors(m, form_data, i)
        ir.append({"AK": AK,
                   "representation": "tensor",
                   "integral_type": "interior_facet_integral",
                   "form_id": form_id,
                   "sub_domain": i,
                   "geometric_dimension": form_data.geometric_dimension,
                   "num_facets": form_data.num_facets})

    return ir

def _extract_tensor_integrals(form, form_data):
    "Extract form containing only tensor representation integrals."
    new_form = Form([])
    for integral in form.integrals():
        if form_data.metadata[integral]["representation"] == "tensor":
            # Get quadrature order and create new integral attaching the order
            # as metadata such that the monomial integration will be aware of
            # quadrature_degree specified by the user on the command line or in forms
            quadrature_degree = form_data.metadata[integral]["quadrature_degree"]
            metadata = {"quadrature_degree": quadrature_degree}
            measure = integral.measure().reconstruct(metadata=metadata)
            integral = Integral(integral.integrand(), measure)
            new_form += Form([integral])
    return new_form

def _compute_cell_tensor(monomial_form, form_data, sub_domain):
    "Compute representation of cell tensor."

    # Extract cell integrals
    monomial_form = _extract_integrals(monomial_form,
                                       form_data,
                                       Measure.CELL,
                                       sub_domain)

    # Compute sum of tensor representations
    terms = _compute_terms(monomial_form, Measure.CELL, None, None)

    return terms

def _compute_exterior_facet_tensors(monomial_form, form_data, sub_domain):
    "Compute representation of exterior facet tensors."

    # Extract exterior facet integrals
    monomial_form = _extract_integrals(monomial_form,
                                       form_data,
                                       Measure.EXTERIOR_FACET,
                                       sub_domain)

    # Compute sum of tensor representations for each facet
    num_facets = form_data.num_facets
    terms = [None for i in range(num_facets)]
    for i in range(num_facets):
        terms[i] = _compute_terms(monomial_form,
                                  Measure.EXTERIOR_FACET,
                                  i, None)

    return terms

def _compute_interior_facet_tensors(monomial_form, form_data, sub_domain):
    "Compute representation of interior facet tensors."

    # Extract interior facet integrals
    monomial_form = _extract_integrals(monomial_form,
                                       form_data,
                                       Measure.INTERIOR_FACET,
                                       sub_domain)

    # Compute sum of tensor representations for each facet-facet pair
    num_facets = form_data.num_facets
    terms = [[None for j in range(num_facets)] for i in range(num_facets)]
    for i in range(num_facets):
        for j in range(num_facets):
            terms[i][j] = _compute_terms(monomial_form,
                                         Measure.INTERIOR_FACET,
                                         i, j)
            reorder_entries(terms[i][j])

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

        # Get quadrature order and pass it on to monomial integration
        quadrature_degree = measure.metadata()["quadrature_degree"]

        # Iterate over monomials of integrand
        for monomial in integrand.monomials:

            # Compute reference tensor
            A0 = ReferenceTensor(monomial,
                                 domain_type,
                                 facet0, facet1,
                                 quadrature_degree)

            # Compute geometry tensor
            GK = GeometryTensor(monomial)

            # Append term
            terms.append((A0, GK))

    return terms
