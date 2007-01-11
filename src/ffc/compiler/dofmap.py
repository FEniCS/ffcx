__author__ = "Robert C. Kirby (kirby@cs.uchicago.edu) and Anders Logg (logg@simula.no)"
__date__ = "2005-05-03 -- 2007-01-19"
__copyright__ = "Copyright (C) 2005-2007 Kirby/Logg"
__license__  = "GNU GPL Version 2"

# This file was previously named dofmap.py

# FIAT modules
from FIAT.dualbasis import *
from FIAT.shapes import *

# FFC modules
from declaration import *
from alignment import *

# FIXME: Should not be DOLFIN-specific
format = { ("entity", 2, 0) : lambda i : "cell.entities(0)[%d]" % i,
           ("entity", 2, 1) : lambda i : "cell.entities(1)[%d]" % i,
           ("entity", 2, 2) : lambda i : "cell.index()",
           ("entity", 2, 3) : lambda i : "not defined",
           ("entity", 3, 0) : lambda i : "cell.entities(0)[%d]" % i,
           ("entity", 3, 1) : lambda i : "cell.entities(1)[%d]" % i,
           ("entity", 3, 2) : lambda i : "cell.entities(2)[%d]" % i,
           ("entity", 3, 3) : lambda i : "cell.index()",
           ("num",    2, 0) : "mesh.topology().size(0)",
           ("num",    2, 1) : "mesh.topology().size(1)",
           ("num",    2, 2) : "mesh.topology().size(2)",
           ("num",    2, 3) : "not defined",
           ("num",    3, 0) : "mesh.topology().size(0)",
           ("num",    3, 1) : "mesh.topology().size(1)",
           ("num",    3, 2) : "mesh.topology().size(2)",
           ("num",    3, 3) : "mesh.topology().size(3)",
           ("check",  2, 0) : lambda i : "not defined",
           ("check",  2, 1) : lambda i : "cell.alignment(1, %d)" % i,
           ("check",  2, 2) : lambda i : "not defined",
           ("check",  2, 3) : lambda i : "not defined",
           ("check",  3, 0) : lambda i : "not defined",
           ("check",  3, 1) : lambda i : "cell.alignment(1, %d)" % i,
           ("check",  3, 2) : lambda i : "cell.alignment(2, %d)" % i,
           ("check",  3, 3) : lambda i : "not defined" }

class DofMap:

    """A NodeMap maps the nodes (degrees of freedom) on a local element
    to global degrees of freedom."""

    def __init__(self, elements):
        "Create NodeMap."

        # FIXME: Argument should be a single element
        self.local_dimension = elements.space_dimension()

        # Make sure we have a list of elements
        if not isinstance(elements, list):
            elements = [elements]

        # Iterate over elements (handles mixed elements)
        self.declarations = []
        data = { "local offset" : 0, "global offset" : None, "num offset" : 0, "num alignment" : 0 }
        for i in range(len(elements)):
            (declarations, data) = self.compute_nodemap(elements[i], data, i)
            self.declarations += declarations

    def write_edge_reordering(self, nodes_per_entity, subelement):
        "Write table for ordering of nodes on edges."
        num_nodes = nodes_per_entity[1]
        map = edge_reordering(num_nodes)
        name = "static unsigned int edge_reordering_%d[2][%d]" % (subelement, num_nodes)
        rows = ["{" + ", ".join(["%d" % node for node in map[i]]) + "}" for i in range(2)]
        value = "{" + ", ".join(rows) + "}"
        return Declaration(name, value)

    def write_face_reordering(self, nodes_per_entity, subelement):
        "Write table for ordering of nodes on faces."
        num_nodes = nodes_per_entity[2]
        map = face_reordering(num_nodes)
        name = "static unsigned int face_reordering_%d[6][%d]" % (subelement, num_nodes)
        rows = ["{" + ", ".join(["%d" % node for node in map[i]]) + "}" for i in range(6)]
        value = "{" + ", ".join(rows) + "}"
        return Declaration(name, value)

    def write_map(self, nodes_per_entity, entity_ids, num_dims, dim, entity, node, data, subelement):
        "Write map from local to global node."
        local_node = entity_ids[dim][entity][node]
        if nodes_per_entity[dim] == 1:
            global_node = format[("entity", num_dims - 1, dim)](entity)
        elif nodes_per_entity[dim] > 1:
            global_node = reorder(nodes_per_entity, num_dims, dim, entity, node, subelement)
        else:
            raise RuntimeError, "Should be at least one node per entity."
        name = "nodes[%d]" % (data["local offset"] + local_node)
        if data["num offset"] == 0:
            value = global_node 
        else:
            value = "offset + " + global_node
        return Declaration(name, value)

    def write_alignment(self, num_dims, dim, entity, data):
        "Write alignment for given entity."
        if data["num alignment"] == 0:
            name = "int alignment"
        else:
            name = "alignment"
        value = format[("check", num_dims - 1, dim)](entity)
        data["num alignment"] += 1
        return Declaration(name, value)

    def reorder(self, nodes_per_entity, num_dims, dim, entity, node, subelement):
        "Compute reordering of map from local to global node."
        start = "%d*%s" % (nodes_per_entity[dim], format[("entity", num_dims - 1, dim)](entity))
        num_nodes = nodes_per_entity[dim]
        # Don't worry about internal nodes
        if dim == (num_dims - 1):
            return start + " + %d" % node
        # Reorder nodes according to orientation of entity
        if dim == 0:
            raise RuntimeError, "Don't know how to reorder nodes for topological dimension 0 (vertices)."
        elif dim == 1:
            return start + " + edge_reordering_%d[alignment][%d]" % (subelement, node)
        elif dim == 2:
            return start + " + face_reordering_%d[alignment][%d]" % (subelement, node)
        else:
            raise RuntimeError, "Don't know how to reorder nodes for topological dimension 3."

    def write_offset(self, data):
        "Write new offset for global nodes."
        increment = data["global offset"]
        if data["num offset"] == 0:
            name = "int offset"
            value = increment
        else:
            name = "offset"
            value = "offset + " + increment
        data["num offset"] += 1
        return Declaration(name, value)

    def compute_offset(self, nodes_per_entity, num_dims, dim):
        "Compute increment for global offset."
        increment = None
        if nodes_per_entity[dim] == 1:
            increment = format[("num", num_dims - 1, dim)]
        elif nodes_per_entity[dim] > 1:
            increment = "%d*%s" % (nodes_per_entity[dim], format[("num", num_dims - 1, dim)])
        return increment

    def compute_nodemap(self, element, data, subelement):
        "Compute nodemap for given element."

        # Get shape and dual basis from FIAT
        shape = element.fiat_shape
        dual_basis = element.fiat_dual
    
        # Get entity IDs
        entity_ids = dual_basis.entity_ids
        
        # Number of topological dimensions
        num_dims = element.shapedim() + 1

        # Count the entities associated with each topological dimension
        num_entities = [len(entity_range(shape, dim)) for dim in range(num_dims)]

        # Count the nodes associated with each entity
        nodes_per_entity = []
        for dim in range(num_dims):
            if entity_ids[dim]:
                nodes_per_entity += [len(entity_ids[dim][0])]
            else:
                nodes_per_entity += [0]

        # Count the total number of nodes (for a scalar element)
        num_nodes = 0
        for dim in range(num_dims):
            num_nodes += num_entities[dim] * nodes_per_entity[dim]

        # Get the number of vector components
        num_components = dual_basis.num_reps

        # Write table for reordering of nodes on edges
        declarations = []
        if nodes_per_entity[1] > 1:
            declarations += [write_edge_reordering(nodes_per_entity, subelement)]

        # Write table for reordering of nodes on faces
        if num_dims > 3:
            if nodes_per_entity[2] > 1:
                declarations += [write_face_reordering(nodes_per_entity, subelement)]

        # Iterate over vector components
        for component in range(num_components):
            # Iterate over topological dimensions
            for dim in range(num_dims):

                # Iterate over entities for current topological dimension
                for entity in range(num_entities[dim]):

                    # Write alignment (if any)
                    if nodes_per_entity[dim] > 1 and dim < (num_dims - 1):
                        declarations += [write_alignment(num_dims, dim, entity, data)]

                    # Iterate over the nodes associated with the current entity
                    for node in range(nodes_per_entity[dim]):

                        # Write offset (if any)
                        if not data["global offset"] == None:
                            declarations += [write_offset(data)]
                            data["global offset"] = None

                        # Write map from local to global node
                        declarations += [self.write_map(nodes_per_entity, entity_ids, num_dims, dim, entity, node, data, subelement)]

                # Add to global offset
                if nodes_per_entity[dim] > 0:
                    data["global offset"] = self.compute_offset(nodes_per_entity, num_dims, dim)

            # Add to local offset (only for vector elements)
            data["local offset"] += num_nodes

            #for declaration in declarations:
            #    print declaration.name + " = " + declaration.value

        return (declarations, data)
