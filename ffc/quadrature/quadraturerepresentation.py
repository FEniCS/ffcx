"Quadrature representation class for UFL"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@gmail.com)"
__date__ = "2009-01-07"
__copyright__ = "Copyright (C) 2009-2010 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Anders Logg, 2009.
# Last changed: 2010-01-28

# UFL modules
from ufl.classes import Form, Integral, SpatialDerivative
from ufl.algorithms import extract_unique_elements, extract_type, extract_elements
#from ufl.algorithms import , extract_unique_elements, extract_type

# FFC modules
from ffc.log import ffc_assert, error
from ffc.fiatinterface import create_element, create_quadrature
from ffc.fiatinterface import map_facet_points

def compute_integral_ir(domain_type, domain_id, integrals, metadata, form_data, form_id):
    "Compute intermediate represention of integral."

    # Initialise representation
    ir = {"representation":       "quadrature",
          "domain_type":          domain_type,
          "domain_id":            domain_id,
          "form_id":              form_id,
          "geometric_dimension":  form_data.geometric_dimension,
          "num_facets":           form_data.num_facets}

    sorted_integrals = _sort_integrals(integrals, metadata, form_data)
    integrals_dict, psi_tables, quad_weights = _tabulate_basis(sorted_integrals, domain_type, form_data.num_facets)

    ir["quadrature_weights"]  = quad_weights
    ir["psi_tables"]          = psi_tables
    ir["integrals"]           = integrals_dict

    return ir

def _tabulate_basis(sorted_integrals, domain_type, num_facets):
    "Tabulate the basisfunctions and derivatives."

    # Initialise return values.
    quadrature_weights = {}
    psi_tables = {}
    integrals = {}

    # Loop the quadrature points and tabulate the basis values.
    for pr, integral in sorted_integrals.iteritems():

        # Extract number of points and the rule.
        # TODO: The rule is currently unused because the fiatinterface does not
        # implement support for other rules than those defined in FIAT_NEW
        num_points_per_axis, rule = pr

        # Get all unique elements in integral.
        elements = extract_unique_elements(integral)

        # Create a list of equivalent FIAT elements (with same ordering of elements).
        fiat_elements = [create_element(e) for e in elements]

        # Get cell and facet domains.
        cell_domain = elements[0].cell().domain()
        facet_domain = elements[0].cell().facet_domain()

        # TODO: These safety check could be removed for speed (I think?)
        ffc_assert(all(cell_domain == e.cell().domain() for e in elements), \
                    "The cell shape of all elements MUST be equal: " + repr(elements))
        ffc_assert(all(facet_domain == e.cell().facet_domain() for e in elements), \
                    "The facet shape of all elements MUST be equal: " + repr(elements))

        # Make quadrature rule and get points and weights.
        # FIXME: Make create_quadrature() take a rule argument.
        if domain_type == "cell":
            (points, weights) = create_quadrature(cell_domain, num_points_per_axis)
        elif domain_type == "exterior_facet" or domain_type == "interior_facet":
            (points, weights) = create_quadrature(facet_domain, num_points_per_axis)
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
        derivatives = set(extract_type(integral, SpatialDerivative))

        # Loop derivatives and extract multiple derivatives.
        for d in list(derivatives):
            num_deriv = len(extract_type(d, SpatialDerivative))

            # TODO: Safety check, SpatialDerivative only has one operand,
            # and there should be only one element?!
            elem = extract_elements(d.operands()[0])
            ffc_assert(len(elem) == 1, "SpatialDerivative has more than one element: " + repr(elem))
            elem = elem[0]
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
            elif domain_type == "exterior_facet" or domain_type == "interior_facet":
                psi_tables[len_weights][elements[i]] = {}
                for facet in range(num_facets):
                    psi_tables[len_weights][elements[i]][facet] =\
                        element.tabulate(deriv_order, map_facet_points(points, facet))
            else:
                error("Unknown domain_type: %s" % domain_type)

    return (integrals, psi_tables, quadrature_weights)

def _sort_integrals(integrals, metadata, form_data):
    """Sort integrals according to the number of quadrature points needed per axis.
    Only consider those integrals defined on the given domain."""

    sorted_integrals = {}
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
        integral_metadata = integral.measure().metadata()
        # Override if specified in integral metadata
        if not integral_metadata is None:
            if "quadrature_degree" in integral_metadata:
                degree = integral_metadata["quadrature_degree"]
            if "quadrature_rule" in integral_metadata:
                rule = integral_metadata["quadrature_rule"]

        # FIXME: This could take place somewhere else?
        # Compute the required number of points for each axis (exact integration).
        num_points_per_axis = (degree + 1 + 1) / 2 # integer division gives 2m - 1 >= q.

        # Create form and add to dictionary according to number of points and rule.
        form = Form([Integral(integral.integrand(), integral.measure().reconstruct(metadata={}))])
        if not (num_points_per_axis, rule) in sorted_integrals:
            sorted_integrals[(num_points_per_axis, rule)] = form
        else:
            sorted_integrals[(num_points_per_axis, rule)] += form
    # Extract integrals form forms.
    for key, val in sorted_integrals.items():
        if len(val.integrals()) != 1:
            error("Only expected one integral over one subdomain: %s" % repr(val))
        sorted_integrals[key] = val.integrals()[0]

    return sorted_integrals
