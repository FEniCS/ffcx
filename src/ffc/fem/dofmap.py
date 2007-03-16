__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-01-24 -- 2007-03-16"
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
        self.__entity_dofs            = entity_dofs
        self.__dof_entities           = self.__compute_dof_entities(entity_dofs)
        self.__signature              = "FFC dof map for " + element.signature()
        self.__local_dimension        = element.space_dimension()
        self.__num_dofs_per_dimension = self.__compute_num_dofs_per_dimension(entity_dofs)

    def entity_dofs(self):
        """Return a dictionary mapping the mesh entities of the
        reference cell to the degrees of freedom associated with the
        entity"""
        return self.__entity_dofs

    def dof_entities(self):
        """Return a list of which entnties are associated with each dof"""
        return self.__dof_entities

    def signature(self):
        "Return a string identifying the finite element"
        return self.__signature

    def local_dimension(self):
        "Return the dimension of the local finite element function space"
        return self.__local_dimension

    def needed_mesh_entities(self):
        """Return a tuple of topological dimensions for the mesh
        entities that are needed to compute the dof map"""
        return self.__needed_mesh_entities

    def num_dofs_per_dimension(self):
        """Return a tuple of the number of dofs associated with each
        topological dimension"""
        return self.__num_dofs_per_dimension

    def __compute_dof_entities(self, entity_dofs):
        """Compute a tuple of the topological entities that correspond
        to each dof:"""
        dof_entities = {}
        for dim in entity_dofs:
            for entity in entity_dofs[dim]:
                for dof in entity_dofs[dim][entity]:
                    dof_entities[dof] = (dim, entity)
        return dof_entities

    def __compute_num_dofs_per_dimension(self, entity_dofs):
        """Compute a tuple of the number of dofs associated with each
        topological dimension"""

        # Count the number of dofs associated with each topological dimension
        num_dofs_per_dimension = [0 for dim in entity_dofs]
        for dim in entity_dofs:
            num_dofs = [len(entity_dofs[dim][entity]) for entity in entity_dofs[dim]]
            # Check that the number of dofs is equal for each entity
            if not num_dofs[1:] == num_dofs[:-1]:
                raise RuntimeError, "The number of dofs must be equal for all entities within a topological dimension."
            # The number of dofs is equal so pick the first
            num_dofs_per_dimension[dim] = num_dofs[0]

        return tuple(num_dofs_per_dimension)
    
    def __repr__(self):
        "Pretty print"
        return self.signature()
