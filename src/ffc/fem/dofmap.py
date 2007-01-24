__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-01-24 -- 2007-01-24"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# FIXME: Temporary fix, do this in mixed element
from mixedelement import *

class DofMap:

    """A DofMap represents a description of the degrees of a freedom
    of a finite element space, from which the mapping from local to
    global degrees of freedom can be computed."""

    def __init__(self, element):
        "Create dof map from given finite element"

        if not isinstance(element, MixedElement):
            self.__entity_dofs = element.entity_dofs()
        self.__signature = "FFC dof map for " + element.signature()

    def entity_dofs(self):
        """Return a dictionary mapping the mesh entities of the
        reference cell to the degrees of freedom associated with
        the entity"""
        return self.__entity_dofs

    def signature(self):
        "Return a string identifying the finite element"
        return self.__signature
