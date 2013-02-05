"Quadrature representation class for UFL"

# Copyright (C) 2009-2010 Kristian B. Oelgaard
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
#
# First added:  2009-01-07
# Last changed: 2010-05-18

import numpy

# UFL modules
from ufl.classes import Form, Integral, Grad
from ufl.algorithms import extract_unique_elements, extract_type, extract_elements

# FFC modules
from ffc.log import ffc_assert, info, error
from ffc.fiatinterface import create_element
from ffc.fiatinterface import map_facet_points, reference_cell_vertices
from ffc.fiatinterface import cellname2num_facets, entities_per_dim
from ffc.quadrature.quadraturetransformer import QuadratureTransformer
from ffc.quadrature.optimisedquadraturetransformer import QuadratureTransformerOpt
from ffc.quadrature_schemes import create_quadrature
from ffc.representationutils import needs_oriented_jacobian

def compute_integral_ir(domain_type,
                        domain_id,
                        integrals,
                        metadata,
                        form_data,
                        form_id,
                        parameters):
    "Compute intermediate represention of integral."

    info("Computing quadrature representation")

    # Get some cell properties
    cell = form_data.cell
    cellname = cell.cellname()
    num_facets = cellname2num_facets[cellname]
    num_vertices = entities_per_dim[cell.topological_dimension()][0]

    # Initialise representation
    ir = {"representation":       "quadrature",
          "domain_type":          domain_type,
          "domain_id":            domain_id,
          "form_id":              form_id,
          "geometric_dimension":  form_data.geometric_dimension,
          "topological_dimension":form_data.topological_dimension,
          "num_facets":           num_facets,
          "num_vertices":         num_vertices,
          "needs_oriented":       needs_oriented_jacobian(form_data),
          "geo_consts":           {}}

    # Sort integrals and tabulate basis.
    sorted_integrals = _sort_integrals(integrals, metadata, form_data)
    integrals_dict, psi_tables, quad_weights = \
        _tabulate_basis(sorted_integrals, domain_type, form_data)

    # Create dimensions of primary indices, needed to reset the argument 'A'
    # given to tabulate_tensor() by the assembler.
    prim_idims = []
    for ufl_element in form_data.argument_elements:
        element = create_element(ufl_element)
        prim_idims.append(element.space_dimension())
    ir["prim_idims"] = prim_idims

    # Create optimise parameters.
    optimise_parameters = {"eliminate zeros":     False,
                           "ignore ones":         False,
                           "remove zero terms":   False,
                           "optimisation":        False,
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

    # Save the optisation parameters.
    ir["optimise_parameters"] = optimise_parameters

    # Create transformer.
    if optimise_parameters["optimisation"]:
        transformer = QuadratureTransformerOpt(psi_tables,
                                               quad_weights,
                                               form_data.geometric_dimension,
                                               form_data.topological_dimension,
                                               form_data.function_replace_map,
                                               optimise_parameters)
    else:
        transformer = QuadratureTransformer(psi_tables,
                                            quad_weights,
                                            form_data.geometric_dimension,
                                            form_data.topological_dimension,
                                            form_data.function_replace_map,
                                            optimise_parameters)

    # Add tables for weights, name_map and basis values.
    ir["quadrature_weights"]  = quad_weights
    ir["name_map"] = transformer.name_map
    ir["unique_tables"] = transformer.unique_tables

    # Transform integrals.
    if domain_type == "cell":
        # Compute transformed integrals.
        info("Transforming cell integral")
        transformer.update_facets(None, None)
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
            transformer.update_facets(None, None)
            transformer.update_vertex(i)
            terms[i] = _transform_integrals(transformer, integrals_dict, domain_type)
        ir["unique_tables"] = transformer.unique_tables
        ir["name_map"] = transformer.name_map
    else:
        error("Unhandled domain type: " + str(domain_type))
    ir["trans_integrals"] = terms

    # Save tables map, to extract table names for optimisation option -O.
    ir["psi_tables_map"] = transformer.psi_tables_map
    ir["additional_includes_set"] = transformer.additional_includes_set

    return ir

def _tabulate_basis(sorted_integrals, domain_type, form_data):
    "Tabulate the basisfunctions and derivatives."

    # MER: Note to newbies: this code assumes that each integral in
    # the dictionary of sorted_integrals that enters here, has a
    # unique number of quadrature points ...

    # Initialise return values.
    quadrature_weights = {}
    psi_tables = {}
    integrals = {}

    # Extract some cell info
    cell = form_data.cell
    gdim = cell.geometric_dimension()
    cellname = cell.cellname()
    facet_cellname = cell.facet_cellname()
    num_facets = cellname2num_facets[cellname]
    num_vertices = entities_per_dim[cell.topological_dimension()][0]

    # Loop the quadrature points and tabulate the basis values.
    for pr, integral in sorted_integrals.iteritems():

        # Extract number of points and the rule.
        # TODO: The rule is currently unused because the fiatinterface does not
        # implement support for other rules than those defined in FIAT_NEW
        degree, rule = pr

        # Get all unique elements in integral.
        elements = [form_data.element_replace_map.get(e,e) for e in extract_unique_elements(integral)]

        # Create a list of equivalent FIAT elements (with same
        # ordering of elements).
        fiat_elements = [create_element(e) for e in elements]

        # Make quadrature rule and get points and weights.
        # FIXME: Make create_quadrature() take a rule argument.
        if domain_type == "cell":
            (points, weights) = create_quadrature(cellname, degree, rule)
        elif domain_type == "exterior_facet" or domain_type == "interior_facet":
            (points, weights) = create_quadrature(facet_cellname, degree, rule)
        elif domain_type == "point":
            (points, weights) = ([()], numpy.array([1.0,])) # Will be fixed
        else:
            error("Unknown integral type: " + str(domain_type))

        # Add points and rules to dictionary.
        len_weights = len(weights) # The TOTAL number of weights/points
        # TODO: This check should not be needed, remove later.
        ffc_assert(len_weights not in quadrature_weights, \
                    "This number of points is already present in the weight table: " + repr(quadrature_weights))
        quadrature_weights[len_weights] = (weights, points)

        # Add the number of points to the psi tables dictionary.
        # TODO: This check should not be needed, remove later.
        ffc_assert(len_weights not in psi_tables, \
                    "This number of points is already present in the psi table: " + repr(psi_tables))
        psi_tables[len_weights] = {}

        # Add the integral with the number of points as a key to the return integrals.
        # TODO: This check should not be needed, remove later.
        ffc_assert(len_weights not in integrals, \
                    "This number of points is already present in the integrals: " + repr(integrals))
        integrals[len_weights] = integral

        # TODO: This is most likely not the best way to get the highest
        # derivative of an element.
        # Initialise dictionary of elements and the number of derivatives.
        num_derivatives = dict([(e, 0) for e in elements])
        # Extract the derivatives from the integral.
        derivatives = set(extract_type(integral, Grad))

        # Loop derivatives and extract multiple derivatives.
        for d in list(derivatives):
            num_deriv = len(extract_type(d, Grad))

            # TODO: Safety check, Grad only has one operand,
            # and there should be only one element?!
            elem = extract_elements(d.operands()[0])
            ffc_assert(len(elem) == 1,
                       "Grad has more than one element: " + repr(elem))
            elem = form_data.element_replace_map.get(elem[0],elem[0])
            # Set the number of derivatives to the highest value
            # encountered so far.
            num_derivatives[elem] = max(num_derivatives[elem], num_deriv)

        # Loop FIAT elements and tabulate basis as usual.
        for i, element in enumerate(fiat_elements):
            # Get order of derivatives.
            deriv_order = num_derivatives[elements[i]]

            # Tabulate for different integral types and insert table into
            # dictionary based on UFL elements.
            if domain_type == "cell":
                psi_tables[len_weights][elements[i]] =\
                {None: element.tabulate(deriv_order, points)}
            elif (domain_type == "exterior_facet"
                  or domain_type == "interior_facet"):
                psi_tables[len_weights][elements[i]] = {}
                for facet in range(num_facets):
                    psi_tables[len_weights][elements[i]][facet] = \
                        element.tabulate(deriv_order,
                                         map_facet_points(points, facet))

            elif domain_type == "point":
                psi_tables[len_weights][elements[i]] = {}
                for vertex in range(num_vertices):
                    points = (reference_cell_vertices(cellname)[vertex],)

                    psi_tables[len_weights][elements[i]][vertex] =\
                        element.tabulate(deriv_order, points)

            else:
                error("Unknown domain_type: %s" % domain_type)

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
            sorted_integrands[(degree, rule)] = integral.integrand()
        else:
            sorted_integrands[(degree, rule)] += integral.integrand()

    # Create integrals from accumulated integrands.
    sorted_integrals = {}
    for key, integrand in sorted_integrands.items():
        sorted_integrals[key] = Integral(integrand, domain_type, domain_id, {}, None)
    return sorted_integrals

def _transform_integrals(transformer, integrals, domain_type):
    "Transform integrals from UFL expression to quadrature representation."
    transformed_integrals = []
    for point, integral in integrals.items():
        transformer.update_points(point)
        terms = transformer.generate_terms(integral.integrand(), domain_type)
        transformed_integrals.append((point, terms, transformer.functions,
                                      {}, transformer.coordinate, transformer.conditionals))
    return transformed_integrals
