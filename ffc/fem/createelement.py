__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl) and Anders Logg (logg@simula.no)"
__date__ = "2009-03-06 -- 2009-08-26"
__copyright__ = "Copyright (C) 2009 Kristian B. Oelgaard and Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# UFL modules
from ufl import FiniteElement as UFLFiniteElement
from ufl import MixedElement as UFLMixedElement
from ufl import ElementRestriction as UFLElementRestriction
from ufl import TensorElement as UFLTensorElement

# FFC common modules
from ffc.common.log import debug, error

# FFC fem modules
from finiteelement import FiniteElement as FFCFiniteElement
from mixedelement import MixedElement as FFCMixedElement
from quadratureelement import QuadratureElement as FFCQuadratureElement

# Cache for computed elements
_cache = {}

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
        if domain and domain != ufl_element.domain():
            error("Domain of restriction is not equal to the domain of the already restricted element. %s %s"\
                   %(str(domain), str(ufl_element.domain())))
        # Get element and domain
        domain = ufl_element.domain()
        ufl_element = ufl_element.element()

    # Create equivalent FFC element
    if isinstance(ufl_element, UFLFiniteElement):
        # Special handling for quadrature elements
        if ufl_element.family() == "Quadrature":
            ffc_element = FFCQuadratureElement(ufl_element.cell().domain(), ufl_element.degree(), domain)
        else:
            ffc_element = FFCFiniteElement(ufl_element.family(), ufl_element.cell().domain(), ufl_element.degree(), domain)
    elif isinstance(ufl_element, UFLMixedElement):
        sub_elements = [create_element(e, domain) for e in ufl_element.sub_elements()]
        ffc_element = FFCMixedElement(sub_elements, domain)
        # FIXME: This is just a temporary hack to 'support' tensor elements
        if isinstance(ufl_element, UFLTensorElement):
            ffc_element._rank = len(ufl_element._shape)
    else:
        error("Unable to create equivalent FIAT element: %s" % str(ufl_element))

    # Add element to cache
    _cache[ufl_element_hash] = ffc_element

    return ffc_element

