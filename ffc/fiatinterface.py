__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl) and Anders Logg (logg@simula.no)"
__date__ = "2009-03-06"
__copyright__ = "Copyright (C) 2009 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Garth N. Wells, 2009.
# Modified by Marie Rognes, 2009.
# Last changed: 2009-12-18

# UFL modules
from ufl import FiniteElement as UFLFiniteElement
from ufl import MixedElement as UFLMixedElement
from ufl import ElementRestriction as UFLElementRestriction
from ufl import TensorElement as UFLTensorElement

# FIAT modules
#from FIAT.shapes import LINE, TRIANGLE, TETRAHEDRON
from FIAT_NEW.reference_element import ufc_simplex
from FIAT_NEW.lagrange import Lagrange
from FIAT_NEW.brezzi_douglas_marini import BrezziDouglasMarini
from FIAT_NEW.discontinuous_lagrange import DiscontinuousLagrange

# FFC modules
from log import debug
from log import error

# FFC fem modules
from quadratureelement import QuadratureElement as FFCQuadratureElement

# Cache for computed elements
_cache = {}

# Mapping from domain to dimension
domain2dim = {"vertex": 0,
              "interval": 1,
              "triangle": 2,
              "tetrahedron": 3}

# Mapping from family name to class
family2class = {"Lagrange": Lagrange,
                "Brezzi-Douglas-Marini": BrezziDouglasMarini,
                "Discontinuous Lagrange": DiscontinuousLagrange}

def create_element(ufl_element):

    # Use element from cache if in cache
    if ufl_element in _cache:
        print "Reusing element from cache!"
        return _cache[ufl_element]

    # FIXME: hack to avoid circular importing
    from mixedelement import MixedElement as FFCMixedElement

    # Create element
    if isinstance(ufl_element, UFLMixedElement):
        element = FFCMixedElement(ufl_element)
    else:
        element = create_fiat_element(ufl_element)

    # Store in cache
    _cache[ufl_element] = element

    return element

def create_fiat_element(ufl_element):
    "Create FIAT element corresponding to given UFL finite element."

    ElementClass = family2class[ufl_element.family()]
    reference_cell = ufc_simplex(domain2dim[ufl_element.cell().domain()])
    element = ElementClass(reference_cell, ufl_element.degree())

    return element
