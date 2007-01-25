__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-01-24 -- 2007-01-25"
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

        # FIXME: Temporary fix for mixed elements
        if not isinstance(element, MixedElement):
            self.__entity_dofs = element.entity_dofs()
            
        self.__signature = "FFC dof map for " + element.signature()

        self.__local_dimension = element.space_dimension()

        self.__global_dimension = self.__compute_global_dimension(element.entity_dofs())

    def entity_dofs(self):
        """Return a dictionary mapping the mesh entities of the
        reference cell to the degrees of freedom associated with
        the entity"""
        return self.__entity_dofs

    def signature(self):
        "Return a string identifying the finite element"
        return self.__signature

    def local_dimension(self):
        "Return the dimension of the local finite element function space"
        return self.__local_dimension

    def global_dimension(self):
        """Compute global dimension as a tuple of the number or mesh
        entities of each topological dimension"""
        return self.__global_dimension

    def __compute_global_dimension(self, entity_dofs):
        """Compute global dimension as a tuple of the number or mesh
        entities of each topological dimension"""

        # Count the number of dofs associated with each topological dimension
        dofs_per_dimension = [0 for dim in entity_dofs]
        for dim in entity_dofs:
            num_dofs = [len(entity_dofs[dim][entity]) for entity in entity_dofs[dim]]
            # Check that the number of dofs is equal for each entity
            if not num_dofs[1:] == num_dofs[:-1]:
                raise RuntimeError, "The number of dofs must be equal for all entities within a topological dimension."
            # The number of dofs is equal so pick the first
            # FIXME: The if-case is a bug fix for a BDM bug in FIAT
            if len(num_dofs) == 0:
                dofs_per_dimension[dim] = 0
            else:
                dofs_per_dimension[dim] = num_dofs[0]

        return tuple(dofs_per_dimension)
