__author__ = "Robert C. Kirby (kirby@cs.uchicago.edu) and Anders Logg (logg@tti-c.org)"
__date__ = "2005-05-03 -- 2005-09-26"
__copyright__ = "Copyright (c) 2005 Kirby/Logg"
__license__  = "GNU GPL Version 2"

# FIAT modules
from FIAT.dualbasis import *
from FIAT.shapes import *

# FFC modules
from declaration import *
from alignment import *

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

def write_edge_reordering(nodes_per_entity):
    "Write table for ordering of dofs on edges."
    num_nodes = nodes_per_entity[1]
    map = edge_reordering(num_nodes)
    name = "static unsigned int edge_reordering[2][%d]" % num_nodes
    rows = ["{" + ", ".join(["%d" % dof for dof in map[i]]) + "}" for i in range(2)]
    value = "{" + ", ".join(rows) + "}"
    return Declaration(name, value)

def write_face_reordering(nodes_per_entity):
    "Write table for ordering of dofs on faces."
    num_nodes = nodes_per_entity[2]
    map = face_reordering(num_nodes)
    name = "static unsigned int face_reordering[6][%d]" % num_nodes
    rows = ["{" + ", ".join(["%d" % dof for dof in map[i]]) + "}" for i in range(6)]
    value = "{" + ", ".join(rows) + "}"
    return Declaration(name, value)

def write_map(nodes_per_entity, entity_ids, num_dims, dim, entity, node, data):
    "Write map from local to global dof."
    local_dof  = entity_ids[dim][entity][node]
    if nodes_per_entity[dim] == 1:
        global_dof = format[("entity", num_dims - 1, dim)](entity)
    elif nodes_per_entity[dim] > 1:
        global_dof = reorder(nodes_per_entity, num_dims, dim, entity, node)
    else:
        raise RuntimeError, "Should be at least one node per entity."
    name = "dofs[%d]" % (data["local offset"] + local_dof)
    if data["num offset"] == 0:
        value = global_dof 
    else:
        value = "offset + " + global_dof
    return Declaration(name, value)

def write_alignment(num_dims, dim, entity, data):
    "Write alignment for given entity."
    if data["num alignment"] == 0:
        name = "int alignment"
    else:
        name = "alignment"
    value = format[("check", num_dims - 1, dim)](entity)
    data["alignment"] += 1
    return Declaration(name, value)

def reorder(nodes_per_entity, num_dims, dim, entity, node):
    "Compute reordering of map from local to global dof."
    start = "%d*%s" % (nodes_per_entity[dim], format[("entity", num_dims - 1, dim)](entity))
    num_nodes = nodes_per_entity[dim]
    # Don't worry about internal dofs
    if dim == (num_dims - 1):
        return start + " + %d" % node
    # Reorder dofs according to orientation of entity
    if dim == 0:
        raise RuntimeError, "Don't know how to reorder dofs for topological dimension 0 (vertices)."
    elif dim == 1:
        return start + " + edge_reordering[alignment][%d]" % node
    elif dim == 2:
        return start + " + face_reordering[alignment][%d]" % node
    else:
        raise RuntimeError, "Don't know how to reorder dofs for topological dimension 3."

def write_offset(data):
    "Write new offset for global dofs."
    increment = data["global offset"]
    if data["num offset"] == 0:
        name = "int offset"
        value = increment
    else:
        name = "offset"
        value = "offset + " + increment
    data["num offset"] += 1
    return Declaration(name, value)

def compute_offset(nodes_per_entity, num_dims, dim):
    "Compute increment for global offset."
    increment = None
    if nodes_per_entity[dim] == 1:
        increment = format[("num", num_dims - 1, dim)]
    elif nodes_per_entity[dim] > 1:
        increment = "%d*%s" % (nodes_per_entity[dim], format[("num", num_dims - 1, dim)])
    return increment

def compute_dofmap(element, data):
    "Compute dofmap for given element."

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

    # Write table for reordering of dofs on edges
    declarations = []
    if nodes_per_entity[1] > 1:
        declarations += [write_edge_reordering(nodes_per_entity)]

    # Write table for reordering of dofs on faces
    if num_dims > 3:
        if nodes_per_entity[2] > 1:
            declarations += [write_face_reordering(nodes_per_entity)]

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

                    # Write map from local to global dof
                    declarations += [write_map(nodes_per_entity, entity_ids, num_dims, dim, entity, node, data)]

            # Add to global offset
            if nodes_per_entity[dim] > 0:
                data["global offset"] = compute_offset(nodes_per_entity, num_dims, dim)

        # Add to local offset (only for vector elements)
        data["local offset"] += num_nodes

    #for declaration in declarations:
    #    print declaration.name + " = " + declaration.value

    return (declarations, data)

class DofMap:

    """A DofMap maps the dofs (degrees of freedom) on a local element
    to global degrees of freedom."""

    def __init__(self, elements):
        "Create DofMap."

        # Make sure we have a list of elements
        if not isinstance(elements, list):
            elements = [elements]

        # Iterate over elements (handles mixed elements)
        self.declarations = []
        data = { "local offset" : 0, "global offset" : None, "num offset" : 0, "num alignment" : 0 }
        for element in elements:
            (declarations, data) = compute_dofmap(element, data)
            self.declarations += declarations
