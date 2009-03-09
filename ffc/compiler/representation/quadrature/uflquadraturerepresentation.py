"Quadrature representation class"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2009-01-07 -- 2009-02-25"
__copyright__ = "Copyright (C) 2009 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# FFC common modules
from ffc.common.debug import *

# FFC language modules
#from ffc.compiler.language.integral import *

# FFC fem modules
from ffc.fem.quadrature import *
#from ffc.fem.finiteelement import FiniteElement as FIATFiniteElement
#from ffc.fem.quadratureelement import QuadratureElement
#from ffc.fem.vectorelement import VectorElement as FIATVectorElement
#from ffc.fem.vectorelement import VectorQuadratureElement
#from ffc.fem.mixedelement import MixedElement as FIATMixedElement
from ffc.fem.referencecell import map_to_facet
#from ffc.fem.quadratureelement import *
from ffc.fem.createelement import *

# FFC quadrature representation modules
#from elementtensor import *

#from factorization import *
#from tensorreordering import *

try:
    from ufl.classes import FiniteElement, MixedElement, VectorElement, TensorElement, FiniteElementBase, Form, Integral
    from ufl.algorithms.analysis import extract_elements, extract_unique_elements, extract_type
#    from ufl.algorithms.ad import expand_derivatives
#    from ufl.algorithms.graph import *
    from ufl.algorithms.printing import tree_format
    from ufl.algorithms.transformations import *
    from ufl.differentiation import SpatialDerivative

    from ufl.integral import Measure
except:
    pass

#def create_fiat_element(ufl_e):

#    if isinstance(ufl_e, TensorElement):
##            print "tensor: "
#        raise RuntimeError(ufl_e, "Tensor element not supported yet.")
#    elif isinstance(ufl_e, VectorElement):
##            print "vector: "
#        # Handle QuadratureElement separately
#        if ufl_e.family() == "Quadrature":
#            num_points_per_axis = (ufl_e.degree() + 1 + 1) / 2 # integer division gives 2m - 1 >= q
#            return VectorQuadratureElement(ufl_e.cell().domain(), num_points_per_axis, len(ufl_e.sub_elements()))
#        return FIATVectorElement(ufl_e.family(), ufl_e.cell().domain(), ufl_e.degree(), len(ufl_e.sub_elements()))
#    elif isinstance(ufl_e, MixedElement):
##            print "mixed: "
#        sub_elems = [create_fiat_element(e) for e in ufl_e.sub_elements()]
#        return FIATMixedElement(sub_elems)
#    elif isinstance(ufl_e, FiniteElement):
##        print "finite: "
#        # Handle QuadratureElement separately
#        if ufl_e.family() == "Quadrature":
#            # UFL Quadrature element set the degree that the element must be able
#            # to integrate where the FFC QuadratureElement sets the number of
#            # points per axis so we need to compute this number
#            num_points_per_axis = (ufl_e.degree() + 1 + 1) / 2 # integer division gives 2m - 1 >= q
##            print "num_points_per_axis: ", num_points_per_axis
#            return QuadratureElement(ufl_e.cell().domain(), num_points_per_axis)
#        # Default return
#        return FIATFiniteElement(ufl_e.family(), ufl_e.cell().domain(), ufl_e.degree())
#    # Element type not supported (yet?).
#    else:
#        raise RuntimeError(ufl_e, "Unable to create equivalent FIAT element.")
#    return

class QuadratureRepresentation:
    """This class initialises some data structures that are used by the
    quadrature code generator.

    Attributes:

        form                     - ??
        cell_integrals           - UFL integrals {subdomain:{num_quad_points: integral,},}
        exterior_facet_integrals - UFL integrals {subdomain:{num_quad_points: integral,},}
        interior_facet_integrals - UFL integrals {subdomain:{num_quad_points: integral,},}

        psi_tables               - tabulated values of basis functions for all
                                   elements of all integrals. Only the unique
                                   elements are tabulated such that efficiency
                                   is maximised.
        quadrature_weights       - same story as the psi_tables.
        fiat_elements_map        - a dictionary, {ufl_element:fiat_element}
    """

    def __init__(self, form_data):
        "Create tensor representation for given form"

        # Save form
        # TODO: Is this still used? should it be?
        self.form = form_data.form

        # Save useful constants
        # TODO: Is this still used? should it be? The fiat_elements_map can be
        # used instead
#        self.geometric_dimension = form_data.geometric_dimension

        # Initialise tables
        self.psi_tables = {Measure.CELL:{},
                           Measure.EXTERIOR_FACET: {},
                           Measure.INTERIOR_FACET: {}}
        self.quadrature_weights = {Measure.CELL:{},
                           Measure.EXTERIOR_FACET: {},
                           Measure.INTERIOR_FACET: {}}
        # TODO: Check if it is faster to just generate new FIAT elements on the fly
        # Is this still used?
#        self.fiat_elements_map = {}

#        print "\nQR, init, form:\n", self.form
#        print "\nQR, init, form.__repr__():\n", self.form.__repr__()
#        print "\nQR, init, tree_format(form):\n", tree_format(self.form)

        # Get relevant integrals of all types
        cell_integrals = self.__extract_integrals(self.form.cell_integrals())
        exterior_facet_integrals = self.__extract_integrals(self.form.exterior_facet_integrals())
        interior_facet_integrals = self.__extract_integrals(self.form.interior_facet_integrals())

        # Tabulate basis values
        self.cell_integrals = self.__tabulate(cell_integrals)
        # FIXME: Tabulate all facet integrals at the same time?
        self.exterior_facet_integrals = self.__tabulate(exterior_facet_integrals)
        self.interior_facet_integrals = self.__tabulate(interior_facet_integrals)

#        print "\nQR, init, psi_tables:\n", self.psi_tables
#        print "\nQR, init, quadrature_weights:\n", self.quadrature_weights

    def __extract_integrals(self, integrals):
        "Extract relevant integrals for the QuadratureGenerator."
        return [i for i in integrals\
            if i.measure().metadata()["ffc_representation"] == "quadrature"]

    def __sort_integrals_quadrature_order(self, integrals):
        "Sort integrals according to the quadrature order"

        sorted_integrals = {}
        # TODO: We might want to take into account that a form like
        # a = f*g*h*v*u*dx(0, quadrature_order=4) + f*v*u*dx(0, quadrature_order=2),
        # although it involves two integrals of different order, will most
        # likely be integrated faster if one does
        # a = (f*g*h + f)*v*u*dx(0, quadrature_order=4)
        # It will of course only work for integrals defined on the same
        # subdomain and representation
        for integral in integrals:
            order = integral.measure().metadata()["quadrature_order"]

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

    def __tabulate(self, unsorted_integrals):
        "Tabulate the basisfunctions and derivatives."

        return_integrals = {}

        # If we don't get any integrals there's nothing to do
        if not unsorted_integrals:
            return None

        # Sort the integrals according number of points needed per axis to
        # integrate the quadrature order exactly
        sorted_integrals = self.__sort_integrals_quadrature_order(unsorted_integrals)

        # The integral type IS the same for ALL integrals
        integral_type = unsorted_integrals[0].measure().domain_type()

        # Loop the quadrature order and tabulate the basis values
        for num_points_per_axis, form in sorted_integrals.items():

#            if len(form.integrals()) != 1:
#                raise RuntimeError(form, "There should be only one integral at this stage")

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
                raise RuntimeError(shape, "The cell shape of all elements MUST be equal")
            if not all(facet_shape == e.facet_shape() for e in fiat_elements):
                raise RuntimeError(facet_shape, "The facet shape of all elements MUST be equal")

            # Make quadrature rule and get points and weights
            if integral_type == Measure.CELL:
                (points, weights) = make_quadrature(shape, num_points_per_axis)
            elif integral_type == Measure.EXTERIOR_FACET or integral_type == Measure.INTERIOR_FACET:
                (points, weights) = make_quadrature(facet_shape, num_points_per_axis)

            # Add rules to dictionary
            len_weights = len(weights) # The TOTAL number of weights/points
            # TODO: This check should not be needed, remove later
            # Figure out what should happen if two different orders require
            # same number of points (e.g., 2nd and 3rd order both require 4 points (2x2))
            if len_weights in self.quadrature_weights[integral_type]:
                raise RuntimeError(len_weights, "This number of points is already present in the table")
            self.quadrature_weights[integral_type][len_weights] = weights

            # Add the number of points to the psi tables dictionary
            # TODO: This check should not be needed, remove later
            # Figure out what should happen if two different orders require
            # same number of points (e.g., 2nd and 3rd order both require 4 points (2x2))
            if len_weights in self.psi_tables[integral_type]:
                raise RuntimeError(len_weights, "This number of points is already present in the table")
            self.psi_tables[integral_type][len_weights] = {}

            # Sort the integrals according to subdomain and add to the return
            # dictionary
            for i in form.integrals():
                subdomain = i.measure().domain_id()
                if subdomain in return_integrals:
                    if len_weights in return_integrals[subdomain]:
                        raise RuntimeError("There should only be one integral for each number of quadrature points on any given subdomain")
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

            # Loop derivatives and extract multiple derivatives
            for d in list(derivatives):
                num_deriv = len(extract_type(d, SpatialDerivative))

                # TODO: Safety check, SpatialDerivative only has one operand,
                # and there should be only one element?!
                elem = extract_elements(d.operands()[0])
                if not len(elem) == 1:
                    raise RuntimeError
                elem = elem[0]
                # Set the number of derivatives to the highest value
                # encountered so far
                num_derivatives[elem] = max(num_derivatives[elem], num_deriv)

            # Loop FIAT elements and tabulate basis as usual
            for i, element in enumerate(fiat_elements):
                # The order in the two lists (fiat_elements and elements)
                # should be the same

                # Update element map
#                if elements[i] not in self.fiat_elements_map:
#                    self.fiat_elements_map[elements[i]] = element

                # Get order of derivatives
                deriv_order = num_derivatives[elements[i]]

                # Tabulate for different integral types and insert table into
                # dictionary based on UFL elements
                # FIXME: Is restriction needed? I only think it is to pick
                # facet0 or facet1 in case of interior facet integrals.
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

    def __debug(self, i, facet0, facet1):
        "Fancy printing of progress"
        if facet0 == facet1 == None:
            debug("Computing quadrature representation for term %d..." % i)
        elif facet1 == None:
            debug("Computing quadrature representation for facet %d, term %d..." % (facet0, i))
        else:
            debug("Computing quadrature representation for facets (%d, %d), term %d..." % (facet0, facet1, i))
