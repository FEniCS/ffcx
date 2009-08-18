__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl) and Anders Logg (logg@simula.no)"
__date__ = "2009-03-06 -- 2009-03-08"
__copyright__ = "Copyright (C) 2009 Kristian B. Oelgaard and Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# UFL modules
from ufl import FiniteElement as UFLFiniteElement
from ufl import MixedElement as UFLMixedElement
from ufl import ElementRestriction as UFLElementRestriction
from ufl.geometry import as_cell

# FFC common modules
from ffc.common.log import debug

# FFC fem modules
from finiteelement import FiniteElement as FFCFiniteElement
from mixedelement import MixedElement as FFCMixedElement
from quadratureelement import QuadratureElement as FFCQuadratureElement

# Cache for computed elements
_cache = {}

def create_element(ufl_element):
    "Create FFC element (wrapper for FIAT element) from UFL element."

    # Check cache
    if ufl_element in _cache:
        debug("Found element in element cache: " + str(ufl_element))
        return _cache[ufl_element]

    # TODO: This is a hack to get started on restriction. For now we can get
    # around with just having the restriction as a property of a FFC element.
    # If it works out, we could change this in UFL?
    domain = None
    ufl_element_hash = ufl_element
    if isinstance(ufl_element, UFLElementRestriction):
        # Get restriction and overwrite element (which is why I save hash)
        domain = as_cell(ufl_element._domain)
        ufl_element = ufl_element._element
        print "\nelement: ", ufl_element
        print "domain: ", domain

    # Create equivalent FFC element
    if isinstance(ufl_element, UFLFiniteElement):
        # Special handling for quadrature elements
        if ufl_element.family() == "Quadrature":
            return create_quadrature_element(ufl_element)
        ffc_element = FFCFiniteElement(ufl_element.family(), ufl_element.cell().domain(), ufl_element.degree())
    elif isinstance(ufl_element, UFLMixedElement):
        sub_elements = [create_element(e) for e in ufl_element.sub_elements()]
        # If the mixed element is restricted, all sub elements must also be restricted
        if domain:
            for e in sub_elements:
                e.domain = domain
        # FIXME: If any of the sub elements are restricted also restrict the mixed element.
        # This is only needed because evaluate_basis etc. doesn't currently support restrictions.
        for e in sub_elements:
            if e.domain:
                domain = e.domain
        ffc_element = FFCMixedElement(sub_elements)
    else:
        raise RuntimeError, ("Unable to create equivalent FIAT element: %s" % str(ufl_element))

    if domain:
        ffc_element.domain = domain

    # Add element to cache
    _cache[ufl_element_hash] = ffc_element

    return ffc_element

def create_quadrature_element(ufl_element):
    "Create FFC quadrature element from UFL quadrature element."
    # Compute the needed number of points to integrate the polynomial degree
    # integer division gives 2*(num_points) - 1 >= polynomial_degree
    num_points_per_axis = (ufl_element.degree() + 1 + 1) / 2
    ffc_element = FFCQuadratureElement(ufl_element.cell().domain(), num_points_per_axis)
    return ffc_element
