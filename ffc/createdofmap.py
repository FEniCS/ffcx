__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2009-03-08"
__copyright__ = "Copyright (C) 2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian B. Oelgaard
# Last changed: 2009-12-09

# FFC modules
from log import debug
from createelement import create_element
from dofmap import DofMap

# Cache for computed dof maps
_cache = {}

def create_dof_map(ufl_element):
    "Create FFC dof map from UFL element."

    # Check cache
    if ufl_element in _cache:
        debug("Found element in dof map cache: " + str(ufl_element))
        return _cache[ufl_element]

    # Create equivalent FFC element and dof map
    ffc_element = create_element(ufl_element)
    ffc_dof_map = DofMap(ffc_element)

    # Add dof map to cache
    _cache[ufl_element] = ffc_dof_map

    return ffc_dof_map
