__author__ = "Robert C. Kirby (kirby@cs.uchicago.edu) and Anders Logg (logg@tti-c.org)"
__date__ = "2005-05-03 -- 2005-05-11"
__copyright__ = "Copyright (c) 2005 Kirby/Logg"
__license__  = "GNU GPL Version 2"

# FIAT modules
from FIAT.dualbasis import *
from FIAT.shapes import *

# FFC modules
from declaration import *

# FIXME: Should not be DOLFIN-specific
format = { ("entity", 2, 0) : lambda i : "cell.nodeID(%d)" % i,
           ("entity", 2, 1) : lambda i : "cell.edgeID(%d)" % i,
           ("entity", 2, 2) : lambda i : "cell.id()",
           ("entity", 2, 3) : lambda i : "not defined",
           ("entity", 3, 0) : lambda i : "cell.nodeID(%d)" % i,
           ("entity", 3, 1) : lambda i : "cell.edgeID(%d)" % i,
           ("entity", 3, 2) : lambda i : "cell.faceID(%d)" % i,
           ("entity", 3, 3) : lambda i : "cell.id()",
           ("num",    2, 0) : "mesh.noNodes()",
           ("num",    2, 1) : "mesh.noEdges()",
           ("num",    2, 2) : "mesh.noCells()",
           ("num",    2, 3) : "not defined",
           ("num",    3, 0) : "mesh.noNodes()",
           ("num",    3, 1) : "mesh.noEdges()",
           ("num",    3, 2) : "mesh.noFaces()",
           ("num",    3, 3) : "mesh.noCells()",
           ("check",  2, 0) : lambda i : "not defined",
           ("check",  2, 1) : lambda i : "cell.edgeAlignment(%d)" % i,
           ("check",  2, 2) : lambda i : "not defined",
           ("check",  2, 3) : lambda i : "not defined",
           ("check",  3, 0) : lambda i : "not defined",
           ("check",  3, 1) : lambda i : "cell.edgeAlignment(%d)" % i,
           ("check",  3, 2) : lambda i : "cell.faceAlignment(%d)" % i,
           ("check",  3, 3) : lambda i : "not defined" }

class DofMap:

    """A DofMap maps the dofs (degrees of freedom) on a local element
    to global degrees of freedom."""

    def __init__(self, shape, dualbasis):
        "Create DofMap."

        # Get entity IDs
        self.entity_ids = dualbasis.entity_ids
        
        # Number of topological dimensions
        self.num_dims = dimension(shape) + 1

        # Count the entities associated with each topological dimension
        self.num_entities = [len(entity_range(shape, dim)) for dim in range(self.num_dims)]
        #print "num entities: " + str(self.num_entities)

        # Count the nodes associated with each entity
        self.nodes_per_entity = []
        for dim in range(self.num_dims):
            if self.entity_ids[dim]:
                self.nodes_per_entity += [len(self.entity_ids[dim][0])]
            else:
                self.nodes_per_entity += [0]
        #print "nodes per entity: " + str(self.nodes_per_entity)

        # Count the total number of nodes (for a scalar element)
        self.num_nodes = 0
        for dim in range(self.num_dims):
            self.num_nodes += self.num_entities[dim] * self.nodes_per_entity[dim]

        # Get the number of vector components
        self.num_components = dualbasis.num_reps

        self.declarations = []
        count = { "offset" : 0, "alignment" : 0 }
        local_offset  = 0    # Total local offset
        global_offset = None # Current increment for global offset

        # Iterate over vector components
        for component in range(self.num_components):
            # Iterate over topological dimensions
            for dim in range(self.num_dims):

                # Iterate over entities for current topological dimension
                for entity in range(self.num_entities[dim]):

                    # Write alignment (if any)
                    if self.nodes_per_entity[dim] > 1 and dim < (self.num_dims - 1):
                        self.declarations += [self.__write_alignment(dim, entity, count)]

                    # Iterate over the nodes associated with the current entity
                    for node in range(self.nodes_per_entity[dim]):

                        # Write offset (if any)
                        if global_offset:
                            self.declarations += [self.__write_offset(global_offset, count)]
                            global_offset = None

                        # Write map from local to global dof
                        self.declarations += [self.__write_map(dim, entity, node, \
                                                               local_offset, count)]

                # Add to global offset
                if self.nodes_per_entity[dim] > 0:
                    global_offset = self.__compute_offset(dim)

            # Add to local offset (only for vector elements)
            local_offset += self.num_nodes

        #for declaration in self.declarations:
        #    print declaration.name + " = " + declaration.value

    def __write_map(self, dim, entity, node, local_offset, count):
        "Write map from local to global dof."
        local_dof  = self.entity_ids[dim][entity][node]
        if self.nodes_per_entity[dim] == 1:
            global_dof = format[("entity", self.num_dims - 1, dim)](entity)
        elif self.nodes_per_entity[dim] > 1:
            global_dof = self.__reorder(dim, entity, node)
        else:
            raise RuntimeError, "Should be at least one node per entity."
        name = "dofs[%d]" % (local_offset + local_dof)
        if count["offset"] == 0:
            value = global_dof 
        else:
            value = "offset + " + global_dof
        return Declaration(name, value)

    def __write_alignment(self, dim, entity, count):
        "Write alignment for given entity."
        if count["alignment"] == 0:
            name = "int alignment"
        else:
            name = "alignment"
        value = format[("check", self.num_dims - 1, dim)](entity)
        count["alignment"] += 1
        return Declaration(name, value)

    def __reorder(self, dim, entity, node):
        "Compute reordering of map from local to global dof."
        start = "%d*%s" % (self.nodes_per_entity[dim], format[("entity", self.num_dims - 1, dim)](entity))
        # Don't worry about internal dofs
        if dim == (self.num_dims - 1):
            return start + " + %d" % node
        # Reorder dofs according to orientation of entity
        if dim == 0:
            raise RuntimeError, "Don't know how to reorder dofs for topological dimension 0 (vertices)."
        elif dim == 1:
            return start + " + ( alignment == 0 ? %d : %d )" % (node, self.nodes_per_entity[dim] - 1 - node)
        elif dim == 2:
            raise RuntimeError, "Reordering of dofs on faces not yet implemented."
        else:
            raise RuntimeError, "Don't know how to reorder dofs for topological dimension 3."

    def __write_offset(self, global_offset, count):
        "Write new offset for global dofs."
        increment = global_offset
        if count["offset"] == 0:
            name = "int offset"
            value = increment
        else:
            name = "offset"
            value = "offset + " + increment
        count["offset"] += 1
        return Declaration(name, value)

    def __compute_offset(self, dim):
        "Compute increment for global offset."
        increment = None
        if self.nodes_per_entity[dim] == 1:
            increment = format[("num", self.num_dims - 1, dim)]
        elif self.nodes_per_entity[dim] > 1:
            increment = "%d*%s" % (self.nodes_per_entity[dim], format[("num", self.num_dims - 1, dim)])
        return increment
