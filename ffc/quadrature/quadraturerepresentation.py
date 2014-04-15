"Quadrature representation class for UFL"

# Copyright (C) 2009-2014 Kristian B. Oelgaard
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
# Modified by Anders Logg 2009, 2014
# Modified by Martin Alnaes 2013
#
# First added:  2009-01-07
# Last changed: 2014-04-02

import numpy
from collections import defaultdict

# UFL modules
from ufl.classes import Form, Integral
from ufl.sorting import sorted_expr_sum

# FFC modules
from ffc.log import ffc_assert, info, error
from ffc.fiatinterface import create_element
from ffc.fiatinterface import cellname_to_num_entities

from ffc.representationutils import initialize_integral_ir
from ffc.quadrature.tabulate_basis import tabulate_basis
from ffc.quadrature.parameters import parse_optimise_parameters

from ffc.quadrature.quadraturetransformer import QuadratureTransformer
from ffc.quadrature.optimisedquadraturetransformer import QuadratureTransformerOpt

def compute_integral_ir(itg_data,
                        form_data,
                        form_id,
                        parameters):
    "Compute intermediate represention of integral."

    info("Computing quadrature representation")

    # Initialise representation
    ir = initialize_integral_ir("quadrature", itg_data, form_data, form_id)

    # Create and save the optisation parameters.
    ir["optimise_parameters"] = parse_optimise_parameters(parameters)

    # Sort integrals into a dict with quadrature degree and rule as key
    sorted_integrals = sort_integrals(itg_data.integrals,
                                      itg_data.metadata["quadrature_degree"],
                                      itg_data.metadata["quadrature_rule"])

    # Tabulate quadrature points and basis function values in these points
    integrals_dict, psi_tables, quadrature_rules = \
        tabulate_basis(sorted_integrals, form_data, itg_data)

    # Save tables for quadrature weights and points
    ir["quadrature_weights"] = quadrature_rules # TODO: Rename this ir entry to quadrature_rules

    # Create dimensions of primary indices, needed to reset the argument 'A'
    # given to tabulate_tensor() by the assembler.
    ir["prim_idims"] = [create_element(ufl_element).space_dimension()
                        for ufl_element in form_data.argument_elements]

    # Create transformer.
    if ir["optimise_parameters"]["optimisation"]:
        QuadratureTransformerClass = QuadratureTransformerOpt
    else:
        QuadratureTransformerClass = QuadratureTransformer

    transformer = QuadratureTransformerClass(psi_tables,
                                             quadrature_rules,
                                             form_data.geometric_dimension,
                                             itg_data.domain.topological_dimension(),
                                             ir["entitytype"],
                                             form_data.function_replace_map,
                                             ir["optimise_parameters"])

    # Transform integrals.
    cellname = itg_data.domain.cell().cellname()
    ir["trans_integrals"] = _transform_integrals_by_type(ir, transformer, integrals_dict,
                                                         itg_data.integral_type, cellname)

    # Save tables populated by transformer
    ir["name_map"] = transformer.name_map
    ir["unique_tables"] = transformer.unique_tables  # Basis values?

    # Save tables map, to extract table names for optimisation option -O.
    ir["psi_tables_map"] = transformer.psi_tables_map
    ir["additional_includes_set"] = transformer.additional_includes_set

    # Insert empty data which will be populated if optimization is turned on
    ir["geo_consts"] = {}

    # Extract element data for psi_tables, needed for runtime quadrature.
    # This is used by integral type custom_integral.
    ir["element_data"] = _extract_element_data(transformer.element_map)

    return ir

def sort_integrals(integrals, default_quadrature_degree, default_quadrature_rule):
    """Sort and accumulate integrals according to the number of quadrature points needed per axis.

    All integrals should be over the same (sub)domain.
    """

    if not integrals:
        return {}

    # Get domain properties from first integral, assuming all are the same
    integral_type  = integrals[0].integral_type()
    subdomain_id    = integrals[0].subdomain_id()
    domain_label = integrals[0].domain().label()
    domain       = integrals[0].domain() # FIXME: Is this safe? Get as input?
    ffc_assert(all(integral_type == itg.integral_type() for itg in integrals),
               "Expecting only integrals of the same type.")
    ffc_assert(all(domain_label == itg.domain().label() for itg in integrals),
               "Expecting only integrals on the same domain.")
    ffc_assert(all(subdomain_id == itg.subdomain_id() for itg in integrals),
               "Expecting only integrals on the same subdomain.")

    sorted_integrands = defaultdict(list)
    for integral in integrals:
        # Override default degree and rule if specified in integral metadata
        integral_metadata = integral.metadata() or {}
        degree = integral_metadata.get("quadrature_degree", default_quadrature_degree)
        rule = integral_metadata.get("quadrature_rule", default_quadrature_rule)
        assert isinstance(degree, int)
        # Add integrand to dictionary according to degree and rule.
        key = (degree, rule)
        sorted_integrands[key].append(integral.integrand())

    # Create integrals from accumulated integrands.
    sorted_integrals = {}
    for key, integrands in sorted_integrands.items():
        # Summing integrands in a canonical ordering defined by UFL
        integrand = sorted_expr_sum(integrands)
        sorted_integrals[key] = Integral(integrand, integral_type, domain, subdomain_id, {}, None)
    return sorted_integrals

def _transform_integrals_by_type(ir, transformer, integrals_dict, integral_type, cellname):
    num_facets = cellname_to_num_entities[cellname][-2]
    num_vertices = cellname_to_num_entities[cellname][0]

    if integral_type == "cell":
        # Compute transformed integrals.
        info("Transforming cell integral")
        transformer.update_cell()
        terms = _transform_integrals(transformer, integrals_dict, integral_type)

    elif integral_type == "exterior_facet":
        # Compute transformed integrals.
        terms = [None]*num_facets
        for i in range(num_facets):
            info("Transforming exterior facet integral %d" % i)
            transformer.update_facets(i, None)
            terms[i] = _transform_integrals(transformer, integrals_dict, integral_type)

    elif integral_type == "interior_facet":
        # Compute transformed integrals.
        terms = [[None]*num_facets for i in range(num_facets)]
        for i in range(num_facets):
            for j in range(num_facets):
                info("Transforming interior facet integral (%d, %d)" % (i, j))
                transformer.update_facets(i, j)
                terms[i][j] = _transform_integrals(transformer, integrals_dict, integral_type)

    elif integral_type == "point":
        # Compute transformed integrals.
        terms = [None]*num_vertices
        for i in range(num_vertices):
            info("Transforming point integral (%d)" % i)
            transformer.update_vertex(i)
            terms[i] = _transform_integrals(transformer, integrals_dict, integral_type)

    elif integral_type == "custom":

        # Compute transformed integrale: same as for cell integrals
        info("Transforming custom integral")
        transformer.update_cell()
        terms = _transform_integrals(transformer, integrals_dict, integral_type)

    else:
        error("Unhandled domain type: " + str(integral_type))
    return terms

def _transform_integrals(transformer, integrals, integral_type):
    "Transform integrals from UFL expression to quadrature representation."
    transformed_integrals = []
    for point, integral in integrals.items():
        transformer.update_points(point)
        terms = transformer.generate_terms(integral.integrand(), integral_type)
        transformed_integrals.append((point, terms, transformer.function_data,
                                      {}, transformer.coordinate, transformer.conditionals))
    return transformed_integrals

def _extract_element_data(element_map):
    "Extract element data for psi_tables"

    # Iterate over map
    element_data = {}
    for elements in element_map.itervalues():
        for ufl_element, counter in elements.iteritems():

            # Create corresponding FIAT element
            fiat_element = create_element(ufl_element)

            # Compute value size
            shape = ufl_element.value_shape()
            value_size = 1 if shape == () else itertools.product(shape)

            # Store data
            element_data[counter] = {"value_size":      value_size,
                                     "local_dimension": fiat_element.space_dimension()}

    return element_data
