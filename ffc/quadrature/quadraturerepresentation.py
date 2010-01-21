"Quadrature representation class for UFL"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@gmail.com)"
__date__ = "2009-01-07"
__copyright__ = "Copyright (C) 2009-2010 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Anders Logg, 2009.
# Last changed: 2010-01-21

# UFL modules
from ufl.classes import Form, Integral, SpatialDerivative
from ufl.algorithms import extract_unique_elements, extract_type, extract_elements
#from ufl.algorithms import , extract_unique_elements, extract_type

# FFC modules
from ffc.log import ffc_assert
from ffc.fiatinterface import create_element, create_quadrature

def compute_integral_ir(form, form_data, form_id, options):
    "Compute intermediate represention of form integrals."

    # Initialise return value
    ir = []

    # Get quadrature integrals of all types.
    cell_integrals = _extract_quadrature_integrals(form.cell_integrals(), form_data)
    ext_facet_integrals = _extract_quadrature_integrals(form.exterior_facet_integrals(), form_data)
    int_facet_integrals = _extract_quadrature_integrals(form.interior_facet_integrals(), form_data)

    # Compute representation of cell tensors.
    for i in range(form_data.num_cell_domains):
        sorted_integrals = _sort_integrals(cell_integrals, i, form_data)
        psi_tables, quad_weights = _tabulate_basis(sorted_integrals, "cell_integral")
        if psi_tables is None and quad_weights is None:
            continue
        ir.append({"integral_type":       "cell_integral",
                   "representation":      "quadrature",
                   "form_id":             form_id,
                   "sub_domain":          i,
                   "quadrature weights":  quad_weights,
                   "psi_tables":          psi_tables,
                   "integrals":           sorted_integrals})

    # Compute representation of exterior facet tensors.
    for i in range(form_data.num_exterior_facet_domains):
        sorted_integrals = _sort_integrals(ext_facet_integrals, i, form_data)
        psi_tables, quad_weights = _tabulate_basis(sorted_integrals, "exterior_facet_integral")
        if psi_tables is None and quad_weights is None:
            continue
        ir.append({"integral_type":       "exterior_facet_integral",
                   "representation":      "quadrature",
                   "form_id":             form_id,
                   "sub_domain":          i,
                   "quadrature weights":  quad_weights,
                   "psi_tables":          psi_tables,
                   "integrals":           sorted_integrals})

    # Compute representation of interior facet tensors.
    for i in range(form_data.num_interior_facet_domains):
        sorted_integrals = _sort_integrals(int_facet_integrals, i, form_data)
        psi_tables, quad_weights = _tabulate_basis(sorted_integrals, "interior_facet_integral")
        if psi_tables is None and quad_weights is None:
            continue
        ir.append({"integral_type":       "interior_facet_integral",
                   "representation":      "quadrature",
                   "form_id":             form_id,
                   "sub_domain":          i,
                   "quadrature weights":  quad_weights,
                   "psi_tables":          psi_tables,
                   "integrals":           sorted_integrals})

    return ir

def _tabulate_basis(sorted_integrals, integral_type):
    "Tabulate the basisfunctions and derivatives."

    # Initialise return values.
    quadrature_weights = {}
    psi_tables = {}

    # If we don't get any integrals there's nothing to do.
    if not sorted_integrals:
        return (None, None)

    # Loop the quadrature points and tabulate the basis values.
    for pr, form in sorted_integrals.iteritems():

        # Extract number of points and the rule.
        # TODO: The rule is currently unused because the fiatinterface does not
        # implement support for other rules than those defined in FIAT_NEW
        num_points_per_axis, rule = pr

        # Get all unique elements in integrals and convert to list.
        elements = set()
        for i in form.integrals():
            elements.update(extract_unique_elements(i))
        elements = list(elements)

        # Create a list of equivalent FIAT elements.
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
        if integral_type == "cell_integral":
            (points, weights) = create_quadrature(cell_domain, num_points_per_axis)
        elif integral_type == "exterior_facet_integral" or integral_type == "interior_facet_integral":
            (points, weights) = make_quadrature(facet_domain, num_points_per_axis, rule)
        else:
            error("Unknown integral type: " + str(integral_type))

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

        # TODO: This is most likely not the best way to get the highest
        # derivative of an element.
        # Initialise dictionary of elements and the number of derivatives.
        num_derivatives = dict([(e, 0) for e in elements])
        # Extract the derivatives from all integrals.
        derivatives = set()
        for i in form.integrals():
            derivatives.update(extract_type(i, SpatialDerivative))

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
            if integral_type == "cell_integral":
                psi_tables[len_weights][elements[i]] =\
                {None: element.tabulate(deriv_order, points)}
            else:# integral_type == "exterior_facet_integral" or integral_type == "interior_facet_integral":
                psi_tables[len_weights][elements[i]] = {}
                for facet in range(element.cell().num_facets()):
                    psi_tables[len_weights][elements[i]][facet] =\
                        element.tabulate(deriv_order, map_to_facet(cell_domain, points, facet))

    return (quadrature_weights, psi_tables)

def _extract_quadrature_integrals(integrals, form_data):
    "Extract relevant integrals (that needs quadrature representation) for the QuadratureGenerator."
    return [i for i in integrals\
            if form_data.metadata[i]["ffc_representation"] == "quadrature"]

def _sort_integrals(integrals, domain_id, form_data):
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
        # Only include integrals on given subdomain.
        if integral.measure().domain_id() != domain_id:
            continue
        order = form_data.metadata[integral]["quadrature_degree"]
        rule  = form_data.metadata[integral]["quadrature_rule"]

        # FIXME: This could take place somewhere else?
        # Compute the required number of points for each axis (exact integration).
        num_points_per_axis = (order + 1 + 1) / 2 # integer division gives 2m - 1 >= q.

        # Create new form and add to dictionary accordingly.
        form = Form([Integral(integral.integrand(), integral.measure().reconstruct(metadata={}))])
        if not (num_points_per_axis, rule) in sorted_integrals:
            sorted_integrals[(num_points_per_axis, rule)] = form
        else:
            sorted_integrals[(num_points_per_axis, rule)] += form

    return sorted_integrals
