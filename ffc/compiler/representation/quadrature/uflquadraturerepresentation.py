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

    from ufl.integral import Integral as UFLIntegral
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

    def __init__(self, form_data, domain_representations, num_quadrature_points):
        "Create tensor representation for given form"

        # Extract form
        form = form_data.form

        # Save form
        self.form = form

        # Save useful constants
        self.geometric_dimension = form_data.geometric_dimension

        # Set number of specified quadrature points. This value can be set from
        # the command line and is valid for ALL forms in the form file. Default
        # value is zero. It might be obsolete given UFL's capability of
        # specifying the number of points as meta-data of the integral in the
        # forms.
        self.num_quad_points = num_quadrature_points

        # Initialise tables
        self.psi_tables = {UFLIntegral.CELL:{},
                           UFLIntegral.EXTERIOR_FACET: {},
                           UFLIntegral.INTERIOR_FACET: {}}
        self.quadrature_weights = {}

        print "\nQR, init, form:\n", form
        print "\nQR, init, form.__repr__():\n", form.__repr__()

        # Get relevant cell integrals
        cell_integrals = [i for i in form.cell_integrals() if\
            domain_representations[(i.domain_type(), i.domain_id())] == "quadrature"]

        # Attach integrals
        self.cell_integrals = cell_integrals
        self.exterior_facet_integrals = []
        self.interior_facet_integrals = []


#        print "\nQR, cell_integrals:\n", cell_integrals

#        integrand = cell_integrals[0].integrand()
#        print integrand.__doc__
#        print "\nQR, cell_integrals[0].integrand()\n", integrand

#        print "\nQR, integrand.free_indices():\n", integrand.free_indices()
#        print "\nQR, integrand.repeated_indices():\n", integrand.repeated_indices()
#        print "\nQR, integrand.index_dimensions():\n", integrand.index_dimensions()

#        operands = integrand.operands()
#        print "\nQR, integrand.operands():\n", operands
#        for op in operands:
#            print "op: ", op
#            print "op.free_indices(): ", op.free_indices()
#            print "op.repeated_indices(): ", op.repeated_indices()
#            print "op.index_dimensions(): ", op.index_dimensions()

#        basisfunctions = extract_basisfunctions(integrand)
#        print "\nQR, basisfunctions:\n", basisfunctions

#        coefficients = extract_coefficients(integrand)
#        print "\nQR, coefficients:\n", coefficients

#        elements = extract_unique_elements(form)
#        print "\nQR, unique elements:\n", elements
#        element = elements.pop()
#        print "\nQR, ufl element:\n", element
#        print "\nQR, ufl element.value shape:\n", element.value_shape()
#        fem = FiniteElement(element.family(), element.cell().domain(), element.degree())
#        print "\nQR, ffc element:\n", fem


#        derivs = sorted(extract_type(integrand, Derivative), cmp=cmp_counted)
#        derivs = extract_type(integrand, SpatialDerivative)
#        print "\nQR, derivs:\n", derivs
#        for d in derivs:
#            print d.operands()[1]
#            print len(d.operands()[1])
#            print d.operands()[1][0]

#        print "QR, form_data.quad_order: ", form_data.quad_order
#        print "QR, count_nodes(integrand): ", count_nodes(integrand)
#        print "QR, extract_duplications(integrand): ", extract_duplications(integrand)

#        print "compound expansion of integrand"
#        print expand_compounds(integrand)

        # FIXME: The quad_order from form_data is very wrong. The number of
        # quadrature points should be determined for each individual integral.
        self.__tabulate(cell_integrals, form_data.quad_order)
        print "\nQR, init, psi_tables:\n", self.psi_tables
        print "\nQR, init, quadrature_weights:\n", self.quadrature_weights
        # Compute representation of cell tensor
#        self.cell_tensor = self.__compute_cell_tensor(form)
        
        # Compute representation of exterior facet tensors
#        self.exterior_facet_tensors = self.__compute_exterior_facet_tensors(form)

        # Compute representation of interior facet tensors
#        self.interior_facet_tensors = self.__compute_interior_facet_tensors(form)

    def __tabulate(self, integrals, order):
        "Tabulate the basisfunctions and derivatives."

        # FIXME: Get polynomial order for each term and integral, not just one
        # value for entire form, take into account derivatives
        num_points = (order + 1 + 1) / 2 # integer division gives 2m - 1 >= q
        if self.num_quad_points:
            num_points = self.num_quad_points
        self.num_quad_points = num_points

        # The integral type is the same for ALL integrals
        integral_type = integrals[0].domain_type()
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
        if integral_type == UFLIntegral.CELL:
            (points, weights) = make_quadrature(shape, num_points)
        elif integral_type == UFLIntegral.EXTERIOR_FACET or integral_type == UFLIntegral.INTERIOR_FACET:
            (points, weights) = make_quadrature(facet_shape, num_points)

#        print "\npoints: ", points
#        print "\nweights: ", weights
        # Add rules to dict
        len_weights = len(weights)
        self.quadrature_weights[integral_type] = {len_weights: weights}

        # Add the number of points to the psi tables dictionary
        self.psi_tables[integral_type] = {len_weights: {}}


        # FIXME: This is most likely not the best way to get the highest
        # derivative of an element
        num_derivatives = {}
        derivs =  set()
        for i in integrals:
            derivs = derivs | extract_type(i, SpatialDerivative)
        for d in list(derivs):
            num_deriv = len(extract_type(d, SpatialDerivative))
            elem = extract_elements(d.operands()[0])
            # FIXME: there should be only one element?!
            if not len(elem) == 1:
                raise RuntimeError
            elem = elem[0]
            if elem in num_derivatives:
                num_derivatives[elem] = max(num_derivatives[elem], num_deriv)
            else:
                num_derivatives[elem] = num_deriv

        for i, element in enumerate(fiat_elements):
            # The order in the two lists should be the same
            deriv_order = 0
            if num_derivatives:
                deriv_order = num_derivatives[elements[i]]
            # Tabulate for different integral types
            if integral_type == Integral.CELL:
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

    def __compute_cell_tensor(self, form):
        "Compute representation of cell tensor"
        debug_begin("Computing cell tensor")

        # Extract monomials
        monomials = self.__extract_monomials(form, Integral.CELL)
        if len(monomials) == 0:
            debug_end()
            return []

        # Compute factorization
        #FIXME: this will mean something else for quadrature
        factorization = self.__compute_factorization(monomials)

        # Compute sum of tensor representations
        tensors = self.__compute_tensors(monomials, factorization, Integral.CELL, None, None)

        debug_end()
        return tensors

    def __compute_exterior_facet_tensors(self, form):
        "Compute representation of exterior facet tensors"
        debug_begin("Computing exterior facet tensors")

        # Extract monomials
        monomials = self.__extract_monomials(form, Integral.EXTERIOR_FACET)
        if len(monomials) == 0:
            debug_end()
            return []

        # Compute factorization
        #FIXME: this will mean something else for quadrature
        factorization = self.__compute_factorization(monomials)

        # Get the number of facets
        num_facets = form.monomials[0].basisfunctions[0].element.num_facets()

        debug("Number of facets to consider: %d" % num_facets)
        
        # Compute sum of tensor representations for each facet
        tensors = [None for i in range(num_facets)]
        for i in range(num_facets):
            tensors[i] = self.__compute_tensors(monomials, factorization, Integral.EXTERIOR_FACET, i, None)

        debug_end()
        return tensors

    def __compute_interior_facet_tensors(self, form):
        "Compute representation of interior facet tensors"
        debug_begin("Computing interior facet tensors")

        # Extract monomials
        monomials = self.__extract_monomials(form, Integral.INTERIOR_FACET)
        if len(monomials) == 0:
            debug_end()
            return []

        # Compute factorization
        #FIXME: this will mean something else for quadrature
        factorization = self.__compute_factorization(monomials)

        # Get the number of facets
        num_facets = form.monomials[0].basisfunctions[0].element.num_facets()

        debug("Number of facets to consider: %d x %d" % (num_facets, num_facets))
        
        # Compute sum of tensor representations for each facet-facet combination
        tensors = [[None for j in range(num_facets)] for i in range(num_facets)]
        for i in range(num_facets):
            for j in range(num_facets):
                tensors[i][j] = self.__compute_tensors(monomials, factorization, Integral.INTERIOR_FACET, i, j)
#                reorder_entries(terms[i][j])

        debug_end()
        return tensors

    def __extract_monomials(self, form, integral_type):
        "Extract monomials"

        # Extract monomials of given type
        monomials = [m for m in form.monomials if m.integral.type == integral_type]
        if len(monomials) > 0:
            debug("Number of terms to consider: %d" % len(monomials))
        else:
            debug("No terms")

        return monomials

    #FIXME: this will mean something else for quadrature, returns None
    def __compute_factorization(self, monomials):
        "Compute factorization"

        factorization = [None for i in range(len(monomials))]

        num_terms = sum([1 for m in factorization if m == None])
        debug("Number of terms to compute: %d" % num_terms)

        return factorization

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
