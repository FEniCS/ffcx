__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-01-24 -- 2007-03-20"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# Modified by Marie E. Rognes (meg@math.uio.no), 2007

# FIXME: Temporary fix, do this in mixed element
from mixedelement import *

class DofMap:

    """A DofMap represents a description of the degrees of a freedom
    of a finite element space, from which the mapping from local to
    global degrees of freedom can be computed."""

    def __init__(self, element):
        "Create dof map from given finite element"

        # Get entity dofs from element
        entity_dofs = element.entity_dofs()

        # Generate dof map data
        self.__signature              = "FFC dof map for " + element.signature()
        self.__local_dimension        = element.space_dimension()
        self.__entity_dofs            = entity_dofs
        self.__dof_entities           = self.__compute_dof_entities(entity_dofs)

    def signature(self):
        "Return a string identifying the dof map"
        return self.__signature

    def local_dimension(self):
        "Return the dimension of the local finite element function space"
        return self.__local_dimension

    def entity_dofs(self):
        """Return a dictionary mapping the mesh entities of the
        reference cell to the degrees of freedom associated with the
        entity"""
        return self.__entity_dofs

    def dof_entities(self):
        """Return a list of which entnties are associated with each dof"""
        return self.__dof_entities

    def __compute_dof_entities(self, entity_dofs):
        "Compute the entities associated with each dof"
        dof_entities = {}
        offset = 0
        for sub_entity_dofs in entity_dofs:
            for dim in sub_entity_dofs:
                for entity in sub_entity_dofs[dim]:
                    for dof in sub_entity_dofs[dim][entity]:
                        dof_entities[offset + dof] = (dim, entity)
            offset = max(dof_entities) + 1
        return dof_entities

    def __repr__(self):
        "Pretty print"
        return self.signature()
