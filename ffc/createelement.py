__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl) and Anders Logg (logg@simula.no)"
__date__ = "2009-03-06"
__copyright__ = "Copyright (C) 2009 Kristian B. Oelgaard and Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Garth N. Wells 2009
# Last changed: 2009-12-16

# UFL modules
from ufl import FiniteElement as UFLFiniteElement
from ufl import MixedElement as UFLMixedElement
from ufl import ElementRestriction as UFLElementRestriction
from ufl import TensorElement as UFLTensorElement

# FIAT modules
from FIAT.shapes import LINE, TRIANGLE, TETRAHEDRON
from FIAT_NEW.lagrange import Lagrange
from FIAT_NEW.reference_element import ufc_simplex

# FFC modules
from log import debug
from log import error

# FFC fem modules
from mixedelement import MixedElement as FFCMixedElement
from quadratureelement import QuadratureElement as FFCQuadratureElement

# Cache for computed elements
_cache = {}

# Dictionaries of basic element data
ufl_domain_dimension = {"vertex": 0, "interval": 1, "triangle": 2, "tetrahedron": 3}
element_string_to_class = {"Lagrange": Lagrange}

def create_fiat_element(ufl_element, domain):
    # Set up FIAT element:
    Element = element_string_to_class[ufl_element.family()]
    reference_cell = ufc_simplex(ufl_domain_dimension[ufl_element.cell().domain()])
    return Element(reference_cell, ufl_element.degree())

def create_element(ufl_element, domain=None):
    "Create FFC element (wrapper for FIAT element) from UFL element."

    # Check cache
    if ufl_element in _cache:
        debug("Found element in element cache: " + str(ufl_element))
        return _cache[ufl_element]

    # Save the element that we'll use as hash
    ufl_element_hash = ufl_element
    if isinstance(ufl_element, UFLElementRestriction):
        # If we already have a domain, make sure that it is equal to the domain
        # of the restricted element
        if domain and domain != ufl_element.domain_restriction():
            error("Domain of restriction is not equal to the domain of the already restricted element. %s %s"\
                   %(str(domain), str(ufl_element.domain())))
        # Get element and domain
        domain = ufl_element.domain()
        ufl_element = ufl_element.element()

    # Create equivalent FFC element
    if isinstance(ufl_element, UFLFiniteElement):
        if ufl_element.family() == "Quadrature":
            ffc_element = FFCQuadratureElement(ufl_element, domain)
        else:
            #ffc_element = FFCFiniteElement(ufl_element, domain)
            ffc_element = create_fiat_element(ufl_element, domain)

    elif isinstance(ufl_element, UFLMixedElement):
        sub_elements = [create_element(e, domain) for e in ufl_element.sub_elements()]
        ffc_element = FFCMixedElement(sub_elements, repr(ufl_element), ufl_element.value_shape(), domain)
        # FIXME: Temporary hack to 'support' tensor elements
        if isinstance(ufl_element, UFLTensorElement):
            ffc_element._rank = len(ufl_element._shape)
    else:
        error("Unable to create equivalent FIAT element: %s" % str(ufl_element))

    # Add element to cache
    _cache[ufl_element_hash] = ffc_element

    return ffc_element

def extract_elements(element):
    elements = []
    if isinstance(element, FFCMixedElement):
        for i in range(element.num_sub_elements()):
            elements += [extract_elements(element.sub_element(i))]
    else:
        elements += [element]
    return elements
