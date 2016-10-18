# -*- coding: utf-8 -*-
"""This module implements the representation of a multilinear form as
a sum of tensor contractions.

The following possible optimizations are currently not implemented but
might be (re-)implemented in a future version of FFC

  1. Factorization of common reference tensors
"""

# Copyright (C) 2007-2014 Anders Logg
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
# Modified by Martin Sandve Alnæs, 2013

# FFC modules
from ffc.log import info, error
from ffc.representationutils import initialize_integral_ir

# FFC tensor representation modules
from ffc.tensor.monomialextraction import extract_monomial_form
from ffc.tensor.monomialtransformation import transform_monomial_form
from ffc.tensor.referencetensor import ReferenceTensor
from ffc.tensor.geometrytensor import GeometryTensor
from ffc.tensor.tensorreordering import reorder_entries


def compute_integral_ir(itg_data,
                        form_data,
                        form_id,
                        element_numbers,
                        classnames,
                        parameters):
    "Compute intermediate represention of integral."

    info("Computing tensor representation")

    # Extract monomial representation
    integrands = [itg.integrand() for itg in itg_data.integrals]
    monomial_form = extract_monomial_form(integrands, form_data.function_replace_map)

    # Transform monomial form to reference element
    transform_monomial_form(monomial_form)

    # Get some integral properties
    integral_type = itg_data.integral_type
    quadrature_degree = itg_data.metadata["quadrature_degree"]
    quadrature_rule = itg_data.metadata["quadrature_rule"]

    # Get some cell properties
    cell = itg_data.domain.ufl_cell()
    num_facets = cell.num_facets()

    # Helper to simplify code below
    compute_terms = lambda i, j: _compute_terms(monomial_form,
                                                i, j,
                                                integral_type,
                                                quadrature_degree,
                                                quadrature_rule,
                                                cell)

    # Compute representation of cell tensor
    if integral_type == "cell":
        # Compute sum of tensor representations
        terms = compute_terms(None, None)

    elif integral_type == "exterior_facet":
        # Compute sum of tensor representations for each facet
        terms = [compute_terms(i, None) for i in range(num_facets)]

    elif integral_type == "interior_facet":
        # Compute sum of tensor representations for each facet-facet pair
        terms = [[compute_terms(i, j) for j in range(num_facets)] for i in range(num_facets)]
        for i in range(num_facets):
            for j in range(num_facets):
                reorder_entries(terms[i][j])

    else:
        error("Unhandled domain type: " + str(integral_type))

    # Initialize representation and store terms
    ir = initialize_integral_ir("tensor", itg_data, form_data, form_id)
    ir["AK"] = terms

    return ir


def _compute_terms(monomial_form,
                   facet0, facet1,
                   integral_type,
                   quadrature_degree,
                   quadrature_rule,
                   cell):
    "Compute list of tensor contraction terms for monomial form."

    # Compute terms
    terms = []
    for integrand in monomial_form:

        # Iterate over monomials of integrand
        for monomial in integrand.monomials:

            # Compute reference tensor
            A0 = ReferenceTensor(monomial,
                                 integral_type,
                                 facet0, facet1,
                                 quadrature_degree,
                                 quadrature_rule,
                                 cell)

            # Compute geometry tensor
            GK = GeometryTensor(monomial)

            # Append term
            terms.append((A0, GK, None))

    return terms
