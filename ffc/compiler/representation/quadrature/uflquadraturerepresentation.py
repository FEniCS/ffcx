"Quadrature representation class for UFL"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2009-01-07 -- 2009-03-18"
__copyright__ = "Copyright (C) 2009 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Anders Logg, 2009.

# FFC common modules
from ffc.common.log import debug, info, error

# FFC fem modules
from ffc.fem.quadrature import make_quadrature
from ffc.fem.referencecell import map_to_facet
from ffc.fem.createelement import create_element

# UFL modules
from ufl.classes import Form, Integral, SpatialDerivative, Measure
from ufl.algorithms import extract_elements, extract_unique_elements, extract_type

class QuadratureRepresentation:
    """This class initialises some data structures that are used by the
    quadrature code generator.

    Attributes:

        form                     - The entire UFL form
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

    def __init__(self, form_data):
        "Create tensor representation for given form"

        info("Extracting information for quadrature representation.")

        # Save form
        self.form = form_data.form

        # Initialise tables
        self.psi_tables = {Measure.CELL:{},
                           Measure.EXTERIOR_FACET: {},
                           Measure.INTERIOR_FACET: {}}
        self.quadrature_weights = {Measure.CELL:{},
                           Measure.EXTERIOR_FACET: {},
                           Measure.INTERIOR_FACET: {}}

        debug("\nQR, init, form:\n" + str(self.form))

        # Get relevant integrals of all types
        cell_integrals = self.__extract_integrals(self.form.cell_integrals(), form_data)
        exterior_facet_integrals = self.__extract_integrals(self.form.exterior_facet_integrals(), form_data)
        interior_facet_integrals = self.__extract_integrals(self.form.interior_facet_integrals(), form_data)

        # Tabulate basis values
        self.cell_integrals = self.__tabulate(cell_integrals, form_data)
        self.exterior_facet_integrals = self.__tabulate(exterior_facet_integrals, form_data)
        self.interior_facet_integrals = self.__tabulate(interior_facet_integrals, form_data)

        debug("\nQR, init, psi_tables:\n" + str(self.psi_tables))
        debug("\nQR, init, quadrature_weights:\n" + str(self.quadrature_weights))

    def __extract_integrals(self, integrals, form_data):
        "Extract relevant integrals for the QuadratureGenerator."
        return [i for i in integrals\
                if form_data.metadata[i]["ffc_representation"] == "quadrature"]

    def __sort_integrals_quadrature_points(self, integrals, form_data):
        "Sort integrals according to the number of quadrature points needed per axis."

        sorted_integrals = {}
        # TODO: We might want to take into account that a form like
        # a = f*g*h*v*u*dx(0, quadrature_order=4) + f*v*u*dx(0, quadrature_order=2),
        # although it involves two integrals of different order, will most
        # likely be integrated faster if one does
        # a = (f*g*h + f)*v*u*dx(0, quadrature_order=4)
        # It will of course only work for integrals defined on the same
        # subdomain and representation
        for integral in integrals:
            order = form_data.metadata[integral]["quadrature_order"]

            # Compute the required number of points for each axis (exact integration)
            num_points_per_axis = (order + 1 + 1) / 2 # integer division gives 2m - 1 >= q

            # FIXME: This could take place somewhere else. In uflcompiler.py
            # we might want to sort according to number of points rather than
            # order?
            if not num_points_per_axis in sorted_integrals:
                sorted_integrals[num_points_per_axis] = Form([Integral(integral.integrand(), integral.measure().reconstruct(metadata={}))])
            else:
                sorted_integrals[num_points_per_axis] += Form([Integral(integral.integrand(), integral.measure().reconstruct(metadata={}))])

        return sorted_integrals

    def __tabulate(self, unsorted_integrals, form_data):
        "Tabulate the basisfunctions and derivatives."

        # Initialise return values
        return_integrals = {}

        # If we don't get any integrals there's nothing to do
        if not unsorted_integrals:
            return None

        # Sort the integrals according number of points needed per axis to
        # integrate the quadrature order exactly
        sorted_integrals = self.__sort_integrals_quadrature_points(unsorted_integrals, form_data)

        # The integral type IS the same for ALL integrals
        integral_type = unsorted_integrals[0].measure().domain_type()

        # Loop the quadrature points and tabulate the basis values
        for num_points_per_axis, form in sorted_integrals.iteritems():

            # Get all unique elements in integrals and convert to list
            elements = set()
            for i in form.integrals():
                elements.update(extract_unique_elements(i))
            elements = list(elements)

            # Create a list of equivalent FIAT elements
            fiat_elements = [create_element(e) for e in elements]

            # Get shape and facet shape.
            shape = fiat_elements[0].cell_shape()
            facet_shape = fiat_elements[0].facet_shape()

            # TODO: These safety check could be removed for speed (I think?)
            if not all(shape == e.cell_shape() for e in fiat_elements):
                print "elements: ", elements
                error("The cell shape of all elements MUST be equal: ", + str(shape))
            if not all(facet_shape == e.facet_shape() for e in fiat_elements):
                print "elements: ", elements
                error("The facet shape of all elements MUST be equal: " + str(facet_shape))

            # Make quadrature rule and get points and weights
            if integral_type == Measure.CELL:
                (points, weights) = make_quadrature(shape, num_points_per_axis)
            elif integral_type == Measure.EXTERIOR_FACET or integral_type == Measure.INTERIOR_FACET:
                (points, weights) = make_quadrature(facet_shape, num_points_per_axis)
            else:
                error("Unknown integral type: " + str(integral_type))

            # Add rules to dictionary
            len_weights = len(weights) # The TOTAL number of weights/points
            # TODO: This check should not be needed, remove later
            if len_weights in self.quadrature_weights[integral_type]:
                print "weights: ", self.quadrature_weights
                error("This number of points is already present in the weight table: " + str(len_weights))
            self.quadrature_weights[integral_type][len_weights] = weights

            # Add the number of points to the psi tables dictionary
            # TODO: This check should not be needed, remove later
            if len_weights in self.psi_tables[integral_type]:
                print "psi tables: ", self.psi_tables
                error("This number of points is already present in the psi table: " + str(len_weights))
            self.psi_tables[integral_type][len_weights] = {}

            # Sort the integrals according to subdomain and add to the return
            # dictionary
            for i in form.integrals():
                subdomain = i.measure().domain_id()
                if subdomain in return_integrals:
                    if len_weights in return_integrals[subdomain]:
                        error("There should only be one integral for any number of quadrature points on any given subdomain")
                    else:
                        return_integrals[subdomain][len_weights] = i
                else:
                    return_integrals[subdomain] = {len_weights: i}

            # TODO: This is most likely not the best way to get the highest
            # derivative of an element
            # Initialise dictionary of elements and the number of derivatives
            num_derivatives = dict([(e, 0) for e in elements])

            # Extract the derivatives from all integrals
            derivatives = set()
            for i in form.integrals():
                derivatives.update(extract_type(i, SpatialDerivative))
            debug("Derivatives: " + str(derivatives))

            # Loop derivatives and extract multiple derivatives
            for d in list(derivatives):
                num_deriv = len(extract_type(d, SpatialDerivative))

                # TODO: Safety check, SpatialDerivative only has one operand,
                # and there should be only one element?!
                elem = extract_elements(d.operands()[0])
                if not len(elem) == 1:
                    error("SpatialDerivative has more than one element: " + str(elem))
                elem = elem[0]
                # Set the number of derivatives to the highest value
                # encountered so far
                num_derivatives[elem] = max(num_derivatives[elem], num_deriv)
            debug("num_derivatives: " + str(num_derivatives))

            # Loop FIAT elements and tabulate basis as usual
            for i, element in enumerate(fiat_elements):
                # The order in the two lists (fiat_elements and elements)
                # should be the same

                # Get order of derivatives
                deriv_order = num_derivatives[elements[i]]

                # Tabulate for different integral types and insert table into
                # dictionary based on UFL elements
                # FIXME: Is restriction needed? I only think it is necessary to
                # pick facet0 or facet1 in case of interior facet integrals.
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

