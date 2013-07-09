"""This module implements the representation of a multilinear form as
a sum of tensor contractions.

The following possible optimizations are currently not implemented but
might be (re-)implemented in a future version of FFC

  1. Factorization of common reference tensors
  2. FErari optimizations
"""

# Copyright (C) 2007-2013 Anders Logg
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Kristian B. Oelgaard, 2010.
# Modified by Martin Alnaes, 2013
#
# First added:  2007-02-05
# Last changed: 2013-02-10

# UFL modules
from ufl.classes import Form, Measure, Integral

# FFC modules
from ffc.log import info, error
from ffc.representationutils import initialize_integral_ir
from ffc.fiatinterface import cellname_to_num_entities

# FFC tensor representation modules
from ffc.tensor.monomialextraction import extract_monomial_form
from ffc.tensor.monomialextraction import MonomialForm
from ffc.tensor.monomialtransformation import transform_monomial_form
from ffc.tensor.referencetensor import ReferenceTensor
from ffc.tensor.geometrytensor import GeometryTensor
from ffc.tensor.tensorreordering import reorder_entries

def compute_integral_ir(itg_data,
                        form_data,
                        form_id,
                        parameters):
    "Compute intermediate represention of integral."

    info("Computing tensor representation")

    # Extract monomial representation
    monomial_form = extract_monomial_form(itg_data.integrals, form_data.function_replace_map)

    # Transform monomial form to reference element
    transform_monomial_form(monomial_form)

    # Get some cell properties
    cell = form_data.cell
    cellname = cell.cellname()
    facet_cellname = cell.facet_cellname()
    num_facets = cellname_to_num_entities[cellname][-2]

    # Initialize representation
    ir = initialize_integral_ir("tensor", itg_data, form_data, form_id)
    ir["rank"] = form_data.rank

    # Compute representation of cell tensor
    quadrature_degree = itg_data.metadata["quadrature_degree"]
    quadrature_rule = itg_data.metadata["quadrature_rule"]
    if itg_data.domain_type == "cell":

        # Compute sum of tensor representations
        ir["AK"] = _compute_terms(monomial_form,
                                  None, None,
                                  itg_data.domain_type,
                                  quadrature_degree,
                                  quadrature_rule,
                                  cellname,
                                  facet_cellname)

    elif itg_data.domain_type == "exterior_facet":

        # Compute sum of tensor representations for each facet
        terms = [None for i in range(num_facets)]
        for i in range(num_facets):
            terms[i] = _compute_terms(monomial_form,
                                      i, None,
                                      itg_data.domain_type,
                                      quadrature_degree,
                                      quadrature_rule,
                                      cellname,
                                      facet_cellname)
        ir["AK"] = terms

    elif itg_data.domain_type == "interior_facet":

        # Compute sum of tensor representations for each facet-facet pair
        terms = [[None for j in range(num_facets)] for i in range(num_facets)]
        for i in range(num_facets):
            for j in range(num_facets):
                terms[i][j] = _compute_terms(monomial_form,
                                             i, j,
                                             itg_data.domain_type,
                                             quadrature_degree,
                                             quadrature_rule,
                                             cellname,
                                             facet_cellname)
                reorder_entries(terms[i][j])
        ir["AK"] = terms

    else:
        error("Unhandled domain type: " + str(itg_data.domain_type))

    return ir

def _compute_terms(monomial_form,
                   facet0, facet1,
                   domain_type,
                   quadrature_degree,
                   quadrature_rule,
                   cellname,
                   facet_cellname):
    "Compute list of tensor contraction terms for monomial form."

    # Compute terms
    terms = []
    for (integrand, measure) in monomial_form:

        # Iterate over monomials of integrand
        for monomial in integrand.monomials:

            # Compute reference tensor
            A0 = ReferenceTensor(monomial,
                                 domain_type,
                                 facet0, facet1,
                                 quadrature_degree,
                                 quadrature_rule,
                                 cellname,
                                 facet_cellname)

            # Compute geometry tensor
            GK = GeometryTensor(monomial)

            # Append term
            terms.append((A0, GK, None))

    return terms
