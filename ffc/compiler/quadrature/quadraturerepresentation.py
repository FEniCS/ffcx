"Quadrature representation class for UFL"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2009-01-07 -- 2009-08-25"
__copyright__ = "Copyright (C) 2009 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Anders Logg, 2009.

# FFC common modules.
from ffc.common.log import debug, info, ffc_assert, error

# FFC fem modules.
from ffc.fem.quadrature import make_quadrature
from ffc.fem.referencecell import map_to_facet
from ffc.fem.createelement import create_element

# UFL modules.
from ufl.classes import Form, Integral, SpatialDerivative, Measure
from ufl.algorithms import extract_elements, extract_unique_elements, extract_type

class QuadratureRepresentation:
    """This class initialises some data structures that are used by the
    quadrature code generator.

    Attributes:

        num_integrals            - total number of integrals

    Attributes added only when num_integrals is nonzero:

        form_data                - UFL form data
        cell_integrals           - UFL integrals {subdomain:{num_quad_points: integral,},}
        exterior_facet_integrals - UFL integrals {subdomain:{num_quad_points: integral,},}
        interior_facet_integrals - UFL integrals {subdomain:{num_quad_points: integral,},}

        psi_tables               - tabulated values of basis functions and their
                                   derivatives up to the highest encountered order for all
                                   elements of all integrals. Only the unique
                                   elements are tabulated such that efficiency
                                   is maximised.
        quadrature_weights       - same story as the psi_tables.
    """

    def __init__(self, form, form_data):
        "Create tensor representation for given form"

        # Save form.
        self.form_data = form_data

        # Initialise tables.
        self.psi_tables = {Measure.CELL:{},
                           Measure.EXTERIOR_FACET: {},
                           Measure.INTERIOR_FACET: {},
                           Measure.SURFACE: {}}
        self.quadrature_weights = {Measure.CELL:{},
                           Measure.EXTERIOR_FACET: {},
                           Measure.INTERIOR_FACET: {},
                           Measure.SURFACE: {}}

        debug("\nQR, init(), form:\n" + str(form))

        # Get relevant integrals of all types.
        cell_integrals = _extract_integrals(form.cell_integrals(), form_data)
        exterior_facet_integrals = _extract_integrals(form.exterior_facet_integrals(), form_data)
        interior_facet_integrals = _extract_integrals(form.interior_facet_integrals(), form_data)

        # Check number of integrals.
        self.num_integrals = len(cell_integrals) + len(exterior_facet_integrals) + len(interior_facet_integrals)
        if self.num_integrals == 0:
            return

        info("Computing quadrature representation.")

        # Tabulate basis values.
        self.cell_integrals = self.__tabulate(_sort_integrals_quadrature_points(cell_integrals, form_data), Measure.CELL)
        self.exterior_facet_integrals = self.__tabulate(_sort_integrals_quadrature_points(exterior_facet_integrals, form_data), Measure.EXTERIOR_FACET)
        self.interior_facet_integrals = self.__tabulate(_sort_integrals_quadrature_points(interior_facet_integrals, form_data), Measure.INTERIOR_FACET)

        debug("\nQR, init(), psi_tables:\n" + str(self.psi_tables))
        debug("\nQR, init(), quadrature_weights:\n" + str(self.quadrature_weights))

    def __tabulate(self, sorted_integrals, integral_type):
        "Tabulate the basisfunctions and derivatives."

        # Initialise return values.
        return_integrals = {}

        # If we don't get any integrals there's nothing to do.
        if not sorted_integrals:
            return None

        # Loop the quadrature points and tabulate the basis values.
        for pr, form in sorted_integrals.iteritems():

            num_points_per_axis, rule = pr

            # Get all unique elements in integrals and convert to list.
            elements = set()
            for i in form.integrals():
                elements.update(extract_unique_elements(i))
            elements = list(elements)

            # Create a list of equivalent FIAT elements.
            fiat_elements = [create_element(e) for e in elements]

            # Get shape and facet shape.
            shape = fiat_elements[0].cell_shape()
            facet_shape = fiat_elements[0].facet_shape()

            # TODO: These safety check could be removed for speed (I think?)
            ffc_assert(all(shape == e.cell_shape() for e in fiat_elements), \
                       "The cell shape of all elements MUST be equal: " + repr(elements))
            ffc_assert(all(facet_shape == e.facet_shape() for e in fiat_elements), \
                       "The facet shape of all elements MUST be equal: " + repr(elements))

            # Make quadrature rule and get points and weights.
            if integral_type == Measure.CELL:
                (points, weights) = make_quadrature(shape, num_points_per_axis, rule)
            elif integral_type == Measure.EXTERIOR_FACET or integral_type == Measure.INTERIOR_FACET:
                (points, weights) = make_quadrature(facet_shape, num_points_per_axis, rule)
            else:
                error("Unknown integral type: " + str(integral_type))

            # Add rules to dictionary.
            len_weights = len(weights) # The TOTAL number of weights/points
            # TODO: This check should not be needed, remove later.
            ffc_assert(len_weights not in self.quadrature_weights[integral_type], \
                       "This number of points is already present in the weight table: " + repr(self.quadrature_weights))
            self.quadrature_weights[integral_type][len_weights] = (weights, points)

            # Add the number of points to the psi tables dictionary.
            # TODO: This check should not be needed, remove later.
            ffc_assert(len_weights not in self.psi_tables[integral_type], \
                       "This number of points is already present in the psi table: " + repr(self.psi_tables))
            self.psi_tables[integral_type][len_weights] = {}

            # Sort the integrals according to subdomain and add to the return
            # dictionary.
            for i in form.integrals():
                subdomain = i.measure().domain_id()
                if subdomain in return_integrals:
                    ffc_assert(len_weights not in return_integrals[subdomain], \
                               "There should only be one integral for any number of quadrature points on any given subdomain.")
                    return_integrals[subdomain][len_weights] = i
                else:
                    return_integrals[subdomain] = {len_weights: i}

            # TODO: This is most likely not the best way to get the highest
            # derivative of an element.
            # Initialise dictionary of elements and the number of derivatives.
            num_derivatives = dict([(e, 0) for e in elements])

            # Extract the derivatives from all integrals.
            derivatives = set()
            for i in form.integrals():
                derivatives.update(extract_type(i, SpatialDerivative))
            debug("Derivatives: " + str(derivatives))

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
            debug("num_derivatives: " + str(num_derivatives))

            # Loop FIAT elements and tabulate basis as usual.
            for i, element in enumerate(fiat_elements):
                # The order in the two lists (fiat_elements and elements)
                # should be the same.

                # Get order of derivatives.
                deriv_order = num_derivatives[elements[i]]

                # Tabulate for different integral types and insert table into
                # dictionary based on UFL elements.
                if integral_type == Measure.CELL:
                    self.psi_tables[integral_type][len_weights]\
                        [elements[i]] = {None: element.tabulate(deriv_order, points)}
                elif integral_type == Measure.EXTERIOR_FACET:
                    self.psi_tables[integral_type][len_weights][elements[i]] = {}
                    for facet in range(element.num_facets()):
                        self.psi_tables[integral_type][len_weights]\
                            [elements[i]][facet] =\
                            element.tabulate(deriv_order, map_to_facet(points, facet))
                elif integral_type == Measure.INTERIOR_FACET:
                    self.psi_tables[integral_type][len_weights][elements[i]] = {}
                    for facet in range(element.num_facets()):
                        self.psi_tables[integral_type][len_weights]\
                            [elements[i]][facet] =\
                            element.tabulate(deriv_order, map_to_facet(points, facet))

        return return_integrals

def _extract_integrals(integrals, form_data):
    "Extract relevant integrals for the QuadratureGenerator."
    return [i for i in integrals\
            if form_data.metadata[i]["ffc_representation"] == "quadrature"]

def _sort_integrals_quadrature_points(integrals, form_data):
    "Sort integrals according to the number of quadrature points needed per axis."

    sorted_integrals = {}
    # TODO: We might want to take into account that a form like
    # a = f*g*h*v*u*dx(0, quadrature_order=4) + f*v*u*dx(0, quadrature_order=2),
    # although it involves two integrals of different order, will most
    # likely be integrated faster if one does
    # a = (f*g*h + f)*v*u*dx(0, quadrature_order=4)
    # It will of course only work for integrals defined on the same
    # subdomain and representation.
    for integral in integrals:
        order = form_data.metadata[integral]["quadrature_order"]
        rule  = form_data.metadata[integral]["quadrature_rule"]

        # Compute the required number of points for each axis (exact integration).
        num_points_per_axis = (order + 1 + 1) / 2 # integer division gives 2m - 1 >= q

        # FIXME: This could take place somewhere else? In uflcompiler.py perhaps?
        if not (num_points_per_axis, rule) in sorted_integrals:
            sorted_integrals[(num_points_per_axis, rule)] = Form([Integral(integral.integrand(), integral.measure().reconstruct(metadata={}))])
        else:
            sorted_integrals[(num_points_per_axis, rule)] += Form([Integral(integral.integrand(), integral.measure().reconstruct(metadata={}))])

    return sorted_integrals
