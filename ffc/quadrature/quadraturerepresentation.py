"Quadrature representation class for UFL"

# Copyright (C) 2009-2013 Kristian B. Oelgaard
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
# Modified by Anders Logg, 2009.
# Modified by Martin Alnaes, 2013
#
# First added:  2009-01-07
# Last changed: 2013-02-10

import numpy

# UFL modules
from ufl.classes import Form, Integral, Grad
from ufl.algorithms import extract_unique_elements, extract_type, extract_elements
from ufl.sorting import sorted_expr_sum

# FFC modules
from ffc.log import ffc_assert, info, error
from ffc.fiatinterface import create_element
from ffc.fiatinterface import map_facet_points, reference_cell_vertices
from ffc.fiatinterface import cellname_to_num_entities
from ffc.quadrature.quadraturetransformer import QuadratureTransformer
from ffc.quadrature.optimisedquadraturetransformer import QuadratureTransformerOpt
from ffc.quadrature_schemes import create_quadrature
from ffc.representationutils import initialize_integral_ir

def compute_integral_ir(itg_data,
                        form_data,
                        form_id,
                        parameters):
    "Compute intermediate represention of integral."

    info("Computing quadrature representation")

    # Initialise representation
    ir = initialize_integral_ir("quadrature", itg_data, form_data, form_id)

    # Sort integrals into a dict with number of integral points as key
    sorted_integrals = _sort_integrals(itg_data.integrals, itg_data.metadata, form_data)

    # Tabulate quadrature points and basis function values in these points
    integrals_dict, psi_tables, quad_weights = \
        _tabulate_basis(sorted_integrals, itg_data.domain_type, form_data)

    # Save tables for quadrature weights and points
    ir["quadrature_weights"]  = quad_weights

    # Create dimensions of primary indices, needed to reset the argument 'A'
    # given to tabulate_tensor() by the assembler.
    ir["prim_idims"] = [create_element(ufl_element).space_dimension()
                        for ufl_element in form_data.argument_elements]

    # Create and save the optisation parameters.
    ir["optimise_parameters"] = _parse_optimise_parameters(parameters)

    # Create transformer.
    if ir["optimise_parameters"]["optimisation"]:
        QuadratureTransformerClass = QuadratureTransformerOpt
    else:
        QuadratureTransformerClass = QuadratureTransformer
    transformer = QuadratureTransformerClass(psi_tables,
                                             quad_weights,
                                             form_data.geometric_dimension,
                                             form_data.topological_dimension,
                                             ir["entitytype"],
                                             form_data.function_replace_map,
                                             ir["optimise_parameters"])

    # Transform integrals.
    ir["trans_integrals"] = _transform_integrals_by_type(ir, transformer, integrals_dict,
                                                         itg_data.domain_type, form_data.cell)

    # Save tables populated by transformer
    ir["name_map"] = transformer.name_map
    ir["unique_tables"] = transformer.unique_tables  # Basis values?

    # Save tables map, to extract table names for optimisation option -O.
    ir["psi_tables_map"] = transformer.psi_tables_map
    ir["additional_includes_set"] = transformer.additional_includes_set

    # Insert empty data which will be populated if optimization is turned on
    ir["geo_consts"] = {}

    return ir

def _parse_optimise_parameters(parameters):
    optimise_parameters = {"eliminate zeros":     False,
                           "optimisation":        False,
                           "ignore ones":         False,
                           "remove zero terms":   False,
                           "ignore zero tables":  False}

    if parameters["optimize"]:
        optimise_parameters["ignore ones"]        = True
        optimise_parameters["remove zero terms"]  = True
        optimise_parameters["ignore zero tables"] = True

        # Do not include this in below if/else clause since we want to be
        # able to switch on this optimisation in addition to the other
        # optimisations.
        if "eliminate_zeros" in parameters:
            optimise_parameters["eliminate zeros"] = True

        if "simplify_expressions" in parameters:
            optimise_parameters["optimisation"] = "simplify_expressions"
        elif "precompute_ip_const" in parameters:
            optimise_parameters["optimisation"] = "precompute_ip_const"
        elif "precompute_basis_const" in parameters:
            optimise_parameters["optimisation"] = "precompute_basis_const"
        # The current default optimisation (for -O) is equal to
        # '-feliminate_zeros -fsimplify_expressions'.
        else:
            # If '-O -feliminate_zeros' was given on the command line, do not
            # simplify expressions
            if not "eliminate_zeros" in parameters:
                optimise_parameters["eliminate zeros"] = True
                optimise_parameters["optimisation"]    = "simplify_expressions"

    return optimise_parameters

def _transform_integrals_by_type(ir, transformer, integrals_dict, domain_type, cell):
    num_facets = cellname_to_num_entities[cell.cellname()][-2]
    num_vertices = cellname_to_num_entities[cell.cellname()][0]

    if domain_type == "cell":
        # Compute transformed integrals.
        info("Transforming cell integral")
        transformer.update_cell()
        terms = _transform_integrals(transformer, integrals_dict, domain_type)

    elif domain_type == "exterior_facet":
        # Compute transformed integrals.
        terms = [None]*num_facets
        for i in range(num_facets):
            info("Transforming exterior facet integral %d" % i)
            transformer.update_facets(i, None)
            terms[i] = _transform_integrals(transformer, integrals_dict, domain_type)

    elif domain_type == "interior_facet":
        # Compute transformed integrals.
        terms = [[None]*num_facets for i in range(num_facets)]
        for i in range(num_facets):
            for j in range(num_facets):
                info("Transforming interior facet integral (%d, %d)" % (i, j))
                transformer.update_facets(i, j)
                terms[i][j] = _transform_integrals(transformer, integrals_dict, domain_type)

    elif domain_type == "point":
        # Compute transformed integrals.
        terms = [None]*num_vertices
        for i in range(num_vertices):
            info("Transforming point integral (%d)" % i)
            transformer.update_vertex(i)
            terms[i] = _transform_integrals(transformer, integrals_dict, domain_type)
    else:
        error("Unhandled domain type: " + str(domain_type))
    return terms

def _create_quadrature_points_and_weights(domain_type, cell, degree, rule):
    if domain_type == "cell":
        (points, weights) = create_quadrature(cell.cellname(), degree, rule)
    elif domain_type == "exterior_facet" or domain_type == "interior_facet":
        (points, weights) = create_quadrature(cell.facet_cellname(), degree, rule)
    elif domain_type == "point":
        (points, weights) = ([()], numpy.array([1.0,])) # TODO: Will be fixed
    else:
        error("Unknown integral type: " + str(domain_type))
    return (points, weights)

def _find_element_derivatives(expr, elements, element_replace_map):
    "Find the highest derivatives of given elements in expression."
    # TODO: This is most likely not the best way to get the highest
    #       derivative of an element, but it works!

    # Initialise dictionary of elements and the number of derivatives.
    # (Note that elements are already mapped through the element_replace_map)
    num_derivatives = dict((e, 0) for e in elements)

    # Extract the derivatives from the integral.
    derivatives = set(extract_type(expr, Grad))

    # Loop derivatives and extract multiple derivatives.
    for d in list(derivatives):
        # After UFL has evaluated derivatives, only one element
        # can be found inside any single Grad expression
        elem, = extract_elements(d.operands()[0])
        elem = element_replace_map[elem]
        # Set the number of derivatives to the highest value encountered so far.
        num_derivatives[elem] = max(num_derivatives[elem], len(extract_type(d, Grad)))
    return num_derivatives

def domain_to_entity_dim(domain_type, cell):
    tdim = cell.topological_dimension()
    if domain_type == "cell":
        entity_dim = tdim
    elif (domain_type == "exterior_facet" or domain_type == "interior_facet"):
        entity_dim = tdim - 1
    elif domain_type == "point":
        entity_dim = 0
    else:
        error("Unknown domain_type: %s" % domain_type)
    return entity_dim

def _map_entity_points(cell, points, entity_dim, entity):
    # Not sure if this is useful anywhere else than in _tabulate_psi_table!
    tdim = cell.topological_dimension()
    if entity_dim == tdim:
        return points
    elif entity_dim == tdim-1:
        return map_facet_points(points, entity)
    elif entity_dim == 0:
        return (reference_cell_vertices(cell.cellname())[entity],)

def _tabulate_psi_table(domain_type, cell, element, deriv_order, points):
    "Tabulate psi table for different integral types."
    # MSA: I attempted to generalize this function, could this way of
    # handling domain types generically extend to other parts of the code?
    entity_dim = domain_to_entity_dim(domain_type, cell)
    num_entities = cellname_to_num_entities[cell.cellname()][entity_dim]
    psi_table = {}
    for entity in range(num_entities):
        entity_points = _map_entity_points(cell, points, entity_dim, entity)
        # TODO: Use 0 as key for cell and we may be able to generalize other places:
        key = None if domain_type == "cell" else entity
        psi_table[key] = element.tabulate(deriv_order, entity_points)
    return psi_table

def _tabulate_entities(domain_type, cell):
    "Tabulate psi table for different integral types."
    # MSA: I attempted to generalize this function, could this way of
    # handling domain types generically extend to other parts of the code?
    entity_dim = domain_to_entity_dim(domain_type, cell)
    num_entities = cellname_to_num_entities[cell.cellname()][entity_dim]
    entities = set()
    for entity in range(num_entities):
        # TODO: Use 0 as key for cell and we may be able to generalize other places:
        key = None if domain_type == "cell" else entity
        entities.add(key)
    return entities

def insert_nested_dict(root, keys, value):
    for k in keys[:-1]:
        d = root.get(k)
        if d is None:
            d = {}
            root[k] = d
        root = d
    root[keys[-1]] = value

from ufl.classes import CellAvg, FacetAvg

def _tabulate_basis(sorted_integrals, domain_type, form_data):
    "Tabulate the basisfunctions and derivatives."

    # MER: Note to newbies: this code assumes that each integral in
    # the dictionary of sorted_integrals that enters here, has a
    # unique number of quadrature points ...

    # Initialise return values.
    quadrature_weights = {}
    psi_tables = {}
    integrals = {}
    avg_elements = { "cell": [], "facet": [] }

    # Loop the quadrature points and tabulate the basis values.
    for pr, integral in sorted_integrals.iteritems():

        # Extract number of points and the rule.
        degree, rule = pr

        # Get all unique elements in integral.
        ufl_elements = [form_data.element_replace_map[e]
                    for e in extract_unique_elements(integral)]

        # Find all CellAvg and FacetAvg in integrals and extract elements
        for avg, AvgType in (("cell", CellAvg), ("facet", FacetAvg)):
            expressions = extract_type(integral, AvgType)
            avg_elements[avg] = [form_data.element_replace_map[e]
                                 for expr in expressions
                                 for e in extract_unique_elements(expr)]

        # Create a list of equivalent FIAT elements (with same
        # ordering of elements).
        fiat_elements = [create_element(e) for e in ufl_elements]

        # Make quadrature rule and get points and weights.
        (points, weights) = _create_quadrature_points_and_weights(domain_type, form_data.cell, degree, rule)

        # The TOTAL number of weights/points
        len_weights = len(weights)

        # Assert that this is unique
        ffc_assert(len_weights not in quadrature_weights, \
                    "This number of points is already present in the weight table: " + repr(quadrature_weights))
        ffc_assert(len_weights not in psi_tables, \
                    "This number of points is already present in the psi table: " + repr(psi_tables))
        ffc_assert(len_weights not in integrals, \
                    "This number of points is already present in the integrals: " + repr(integrals))

        # Add points and rules to dictionary.
        quadrature_weights[len_weights] = (weights, points)

        # Add the number of points to the psi tables dictionary.
        psi_tables[len_weights] = {}

        # Add the integral with the number of points as a key to the return integrals.
        integrals[len_weights] = integral

        # Find the highest number of derivatives needed for each element
        num_derivatives = _find_element_derivatives(integral.integrand(), ufl_elements,
                                                    form_data.element_replace_map)

        # Loop FIAT elements and tabulate basis as usual.
        for i, element in enumerate(fiat_elements):
            # Tabulate table of basis functions and derivatives in points
            psi_table = _tabulate_psi_table(domain_type, form_data.cell, element,
                                        num_derivatives[ufl_elements[i]], points)

            # Insert table into dictionary based on UFL elements. (None=not averaged)
            psi_tables[len_weights][ufl_elements[i]] = { None: psi_table }

    # Loop over elements found in CellAvg and tabulate basis averages
    len_weights = 1
    for avg in ("cell", "facet"):
        # Doesn't matter if it's exterior or interior
        if avg == "cell":
            avg_domain_type = "cell"
        elif avg == "facet":
            avg_domain_type = "exterior_facet"

        for element in avg_elements[avg]:
            fiat_element = create_element(element)

            # Make quadrature rule and get points and weights.
            (points, weights) = _create_quadrature_points_and_weights(avg_domain_type, form_data.cell, element.degree(), "default")
            wsum = sum(weights)

            # Tabulate table of basis functions and derivatives in points
            entity_psi_tables = _tabulate_psi_table(avg_domain_type, form_data.cell, fiat_element, 0, points)
            rank = len(element.value_shape())

            # Hack, duplicating table with per-cell values for each facet in the case of cell_avg(f) in a facet integral
            actual_entities = _tabulate_entities(domain_type, form_data.cell)
            if len(actual_entities) > len(entity_psi_tables):
                assert len(entity_psi_tables) == 1
                assert avg_domain_type == "cell"
                assert "facet" in domain_type
                v, = entity_psi_tables.values()
                entity_psi_tables = dict((e, v) for e in actual_entities)

            for entity, deriv_table in entity_psi_tables.items():
                deriv, = list(deriv_table.keys()) # Not expecting derivatives of averages
                psi_table = deriv_table[deriv]

                if rank:
                    # Compute numeric integral
                    num_dofs, num_components, num_points = psi_table.shape
                    ffc_assert(num_points == len(weights), "Weights and table shape does not match.")
                    avg_psi_table = numpy.asarray([[[numpy.dot(psi_table[j,k,:], weights) / wsum]
                                                   for k in range(num_components)]
                                                   for j in range(num_dofs)])
                else:
                    # Compute numeric integral
                    num_dofs, num_points = psi_table.shape
                    ffc_assert(num_points == len(weights), "Weights and table shape does not match.")
                    avg_psi_table = numpy.asarray([[numpy.dot(psi_table[j,:], weights) / wsum] for j in range(num_dofs)])

                # Insert table into dictionary based on UFL elements.
                insert_nested_dict(psi_tables, (len_weights, element, avg, entity, deriv), avg_psi_table)

    return (integrals, psi_tables, quadrature_weights)

def _sort_integrals(integrals, metadata, form_data):
    """Sort integrals according to the number of quadrature points needed per axis.
    Only consider those integrals defined on the given domain."""

    # Get domain type and id
    if integrals:
        domain_type = integrals[0].domain_type()
        domain_id = integrals[0].domain_id()
        ffc_assert(all(domain_type == itg.domain_type() for itg in integrals),
                   "Expecting only integrals on the same subdomain.")
        ffc_assert(all(domain_id == itg.domain_id() for itg in integrals),
                   "Expecting only integrals on the same subdomain.")

    sorted_integrands = {}
    # TODO: We might want to take into account that a form like
    # a = f*g*h*v*u*dx(0, quadrature_order=4) + f*v*u*dx(0, quadrature_order=2),
    # although it involves two integrals of different order, will most
    # likely be integrated faster if one does
    # a = (f*g*h + f)*v*u*dx(0, quadrature_order=4)
    # It will of course only work for integrals defined on the same
    # subdomain and representation.
    for integral in integrals:
        # Get default degree and rule.
        degree = metadata["quadrature_degree"]
        rule  = metadata["quadrature_rule"]

        # Override if specified in integral metadata
        integral_compiler_data = integral.compiler_data()
        if not integral_compiler_data is None:
            if "quadrature_degree" in integral_compiler_data:
                degree = integral_compiler_data["quadrature_degree"]
            if "quadrature_rule" in integral_compiler_data:
                rule = integral_compiler_data["quadrature_rule"]

        # Add integrand to dictionary according to degree and rule.
        if not (degree, rule) in sorted_integrands:
            sorted_integrands[(degree, rule)] = [integral.integrand()]
        else:
            sorted_integrands[(degree, rule)] += [integral.integrand()]

    # Create integrals from accumulated integrands.
    sorted_integrals = {}
    for key, integrands in sorted_integrands.items():
        # Summing integrands in a canonical ordering defined by UFL
        integrand = sorted_expr_sum(integrands)
        sorted_integrals[key] = Integral(integrand, domain_type, domain_id, {}, None)
    return sorted_integrals

def _transform_integrals(transformer, integrals, domain_type):
    "Transform integrals from UFL expression to quadrature representation."
    transformed_integrals = []
    for point, integral in integrals.items():
        transformer.update_points(point)
        terms = transformer.generate_terms(integral.integrand(), domain_type)
        transformed_integrals.append((point, terms, transformer.function_data,
                                      {}, transformer.coordinate, transformer.conditionals))
    return transformed_integrals
