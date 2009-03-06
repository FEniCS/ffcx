__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl) and Anders Logg (logg@simula.no)"
__date__ = "2009-03-06 -- 2009-03-06"
__copyright__ = "Copyright (C) 2009 Kristian B. Oelgaard and Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# UFL modules
from ufl import FiniteElement as UFLFiniteElement
from ufl import MixedElement as UFLMixedElement

# FFC common modules
from ffc.common.log import debug

# FFC fem modules
from finiteelement import FiniteElement as FFCFiniteElement
from mixedelement import MixedElement as FFCMixedElement
from quadratureelement import QuadratureElement as FFCQuadratureElement

# Cache for computed elements
_element_cache = {}

def create_element(ufl_element):
    "Create FFC element (wrapper for FIAT element) from UFL element."

    # Check cache
    if ufl_element in _element_cache:
        debug("Found element in cache:", ufl_element)
        return _element_cache[ufl_element]

    # Special handling for quadrature elements
    if ufl_element.family() == "Quadrature":
        return create_quadrature_element(ufl_element)

    # Create equivalent FFC element
    if isinstance(ufl_element, UFLFiniteElement):
        ffc_element = FFCFiniteElement(ufl_element.family(), ufl_element.cell().domain(), ufl_element.degree())
    elif isinstance(ufl_element, UFLMixedElement):
        sub_elements = [create_element(e) for e in ufl_element.sub_elements()]
        ffc_element = FFCMixedElement(sub_elements)
    else:
        raise RuntimeError, ("Unable to create equivalent FIAT element: %s" % str(ufl_element))

    # Add element to cache
    _element_cache[ufl_element] = ffc_element

    return ffc_element

def create_quadrature_element(ufl_element):
    "Create FFC quadrature element from UFL quadrature element."
    raise NotImplementedError
