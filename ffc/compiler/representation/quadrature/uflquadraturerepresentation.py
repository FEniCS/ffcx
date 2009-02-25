"Quadrature representation class"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2009-01-07 -- 2009-01-08"
__copyright__ = "Copyright (C) 2009 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# FFC common modules
from ffc.common.debug import *

# FFC language modules
#from ffc.compiler.language.integral import *

# FFC fem modules
from ffc.fem.quadrature import *
from ffc.fem.finiteelement import FiniteElement as FIATFiniteElement
from ffc.fem.vectorelement import VectorElement as FIATVectorElement
from ffc.fem.mixedelement import MixedElement as FIATMixedElement
from ffc.fem.referencecell import *
#from ffc.fem.quadratureelement import *

# FFC quadrature representation modules
#from elementtensor import *

#from factorization import *
#from tensorreordering import *

try:
    from ufl.classes import FiniteElement, MixedElement, VectorElement, FiniteElementBase
    from ufl.algorithms.analysis import *
    from ufl.algorithms.graph import *
    from ufl.algorithms.transformations import *
    from ufl.differentiation import SpatialDerivative

    from ufl.integral import Measure
except:
    pass

class QuadratureRepresentation:
    """This class uses quadrature to represent a given multilinear form.

    Attributes:

        form                            - the form generating the quadrature representation
        cell_tensor                     - the representation of the cell tensor
        exterior_facet_tensors          - the representation of the interior facet tensors,
                                          one for each facet
        interior_facet_tensors          - the representation of the exterior facet tensors,
                                          one for each facet-facet combination
        num_user_specified_quad_points  - the number of desired quadrature points specified
                                          by the user. Will be used for ALL terms

    """

    def __init__(self, form_data):
        "Create tensor representation for given form"

        # Extract form
        form = form_data.form

        # Save form
        self.form = form

        # Save useful constants
        self.geometric_dimension = form_data.geometric_dimension

        # Initialise tables
        self.psi_tables = {Measure.CELL:{},
                           Measure.EXTERIOR_FACET: {},
                           Measure.INTERIOR_FACET: {}}
        self.quadrature_weights = {}

        print "\nQR, init, form:\n", form
        print "\nQR, init, form.__repr__():\n", form.__repr__()

        # Get relevant integrals of all types and attach
        self.cell_integrals = self.__extract_integrals(form.cell_integrals())
        self.exterior_facet_integrals =\
                    self.__extract_integrals(form.exterior_facet_integrals())
        self.interior_facet_integrals =\
                    self.__extract_integrals(form.interior_facet_integrals())

        # Tabulate basis values
        print "cell_integrals: ", self.cell_integrals
        self.__tabulate(self.cell_integrals)
#        self.__tabulate(self.exterior_facet_integrals)
#        self.__tabulate(self.interior_facet_integrals)

        print "\nQR, init, psi_tables:\n", self.psi_tables
        print "\nQR, init, quadrature_weights:\n", self.quadrature_weights
        # Compute representation of cell tensor

    def __extract_integrals(self, integrals):
        "Extract relevant integrals for the QuadratureGenerator."
        print "I: ", integrals
        return [i for i in integrals if\
               i.measure().metadata()["ffc"]["representation"] == "quadrature"]

    def __tabulate(self, integrals):
        "Tabulate the basisfunctions and derivatives."

        # FIXME: Get polynomial order for each term and integral, not just one
        # value for entire form, take into account derivatives
        order = integrals[0].measure().metadata()["quadrature_order"]

        num_points = (order + 1 + 1) / 2 # integer division gives 2m - 1 >= q

        # The integral type is the same for ALL integrals
        integral_type = integrals[0].measure().domain_type()
        print "\nQR, tabulate, integral_type:\n", integral_type

        # Get all unique elements in integrals and convert to list
        elements = set()
        for i in integrals:
            elements = elements | extract_unique_elements(i)
        elements = list(elements)
        print "\nQR, tabulate, unique elements:\n", elements
        fiat_elements = [self.__create_fiat_elements(e) for e in elements]
        print "\nQR, tabulate, unique fiat elements:\n", fiat_elements
        print "\nQR, tabulate, unique elements[0].value_shape():\n", elements[0].value_shape()
        print "\nQR, tabulate, unique fiat_elements[0].rank():\n", fiat_elements[0].value_rank()


        # Get shape (check first one, all should be the same)
        # FIXME: Does this still hold true? At least implement check
        shape = fiat_elements[0].cell_shape()
        facet_shape = fiat_elements[0].facet_shape()
#        print "shape: ", shape
#        print "facet_shape: ", facet_shape
        # Make quadrature rule
        # FIXME: need to make different rules for different terms, depending
        # on order
        # Create quadrature rule and get points and weights
        if integral_type == Measure.CELL:
            (points, weights) = make_quadrature(shape, num_points)
        elif integral_type == Measure.EXTERIOR_FACET or integral_type == Measure.INTERIOR_FACET:
            (points, weights) = make_quadrature(facet_shape, num_points)

#        print "\npoints: ", points
#        print "\nweights: ", weights
        # Add rules to dict
        len_weights = len(weights)
        self.quadrature_weights[integral_type] = {len_weights: weights}

        # Add the number of points to the psi tables dictionary
        self.psi_tables[integral_type] = {len_weights: {}}

        # TODO: This is most likely not the best way to get the highest
        # derivative of an element
        # Initialise dictionary of elements and the number of derivatives
        num_derivatives = dict([(e,0) for e in elements])
        derivs =  set()
        for i in integrals:
            derivs = derivs | extract_type(i, SpatialDerivative)
        for d in list(derivs):
            num_deriv = len(extract_type(d, SpatialDerivative))
            elem = extract_elements(d.operands()[0])
            # TODO: there should be only one element?!
            if not len(elem) == 1:
                raise RuntimeError
            elem = elem[0]
            num_derivatives[elem] = max(num_derivatives[elem], num_deriv)

        for i, element in enumerate(fiat_elements):
            # The order in the two lists should be the same
            deriv_order = 0
            if num_derivatives:
                deriv_order = num_derivatives[elements[i]]
            # Tabulate for different integral types
            if integral_type == Measure.CELL:
#                self.psi_tables[integral_type][len_weights][(element, None)] =\
#                                          element.tabulate(deriv_order, points)
                # Create dictionary based on UFL elements
                self.psi_tables[integral_type][len_weights][(elements[i], None)] =\
                                          element.tabulate(deriv_order, points)
        return

    def __create_fiat_elements(self, ufl_e):

        if isinstance(ufl_e, VectorElement):
            return FIATVectorElement(ufl_e.family(), ufl_e.cell().domain(), ufl_e.degree(), len(ufl_e.sub_elements()))
        elif isinstance(ufl_e, MixedElement):
            sub_elems = [self.__create_fiat_elements(e) for e in ufl_e.sub_elements()]
            return FIATMixedElement(sub_elems)
        elif isinstance(ufl_e, FiniteElement):
            return FIATFiniteElement(ufl_e.family(), ufl_e.cell().domain(), ufl_e.degree())
        # Element type not supported (yet?) TensorElement will trigger this.
        else:
            raise RuntimeError(ufl_e, "Unable to create equivalent FIAT element.")
        return


    def __compute_tensors(self, monomials, factorization, integral_type, facet0, facet1):
        "Compute terms and factorize common reference tensors"

        # Compute terms
        num_tensors = len(monomials)
        tensors = [None for i in range(num_tensors)]

        for i in range(num_tensors):

            # Get monomial
            m = monomials[i]

            # Only consider monomials of given integral type
            if not m.integral.type == integral_type:
                continue

            # Compute element tensor
            self.__debug(i, facet0, facet1)
            tensors[i] = ElementTensor(m, facet0, facet1, self.num_user_specified_quad_points)
            debug("done")

        return tensors

    def __debug(self, i, facet0, facet1):
        "Fancy printing of progress"
        if facet0 == facet1 == None:
            debug("Computing quadrature representation for term %d..." % i)
        elif facet1 == None:
            debug("Computing quadrature representation for facet %d, term %d..." % (facet0, i))
        else:
            debug("Computing quadrature representation for facets (%d, %d), term %d..." % (facet0, facet1, i))
