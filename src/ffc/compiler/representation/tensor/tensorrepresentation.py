__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-02-05 -- 2007-02-05"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC common modules
from ffc.common.debug import *

class TensorRepresentation:
    """This class represents a given multilinear form as a tensor
    contraction, or more precisely, a sum of tensor contractions for
    each type of integral: cell, exterior facet and interior facet."""

    def __init__(self, form):
        "Create tensor representation for given form"

        # Compute representation of element tensor
        self.__element_tensor = self.__compute_element_tensor(form)
        
        # Compute representation of exterior facet tensors
        self.__exterior_facet_tensors = self.__compute_exterior_facet_tensors(form)

        # Compute representation of interior facet tensors
        self.__interior_facet_tensors = self.__compute_interior_facet_tensors(form)
        
    def element_tensor():
        "Return representation of element tensor"
        return self.__element_tensor

    def exterior_facet_tensor(self, i, j):
        "Return representation of exterior facet tensor on local facet i"
        return self.__exterior_facet_tensors[i]

    def interior_facet_tensor(self, i):
        """Return representation of interior facet tensor on intersecton
        of local facets i and j"""
        
        return self.__interior_facet_tensors[i][j]

    def __compute_element_tensor(self, form):
        "Compute representation of element tensor"
        debug("Computing element tensor...")
        debug("not implemented")
        return "Not implemented"

    def __compute_exterior_facet_tensors(self, form):
        "Compute representation of exterior facet tensors"
        debug("Computing exterior facet tensors...")
        debug("not implemented")
        return "Not implemented"

    def __compute_interior_facet_tensors(sefl, form):
        "Compute representation of interior facet tensors"
        debug("Computing interior facet tensors...")
        debug("not implemented")
        return "Not implemented"
